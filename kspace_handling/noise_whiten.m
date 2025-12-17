function [kspace_whitened, info] = noise_whiten(kspace, obj, opts) %#ok<INUSD>
%NOISE_WHITEN  Coil noise whitening with robust noise covariance estimation.
%
% Drop-in compatible with noise_whiten(kspace,obj).
% Assumes receiver channels are in the LAST dimension: [Nx Ny ... Nc]
%
% Optional opts fields:
%   opts.method        = 'outer' (default) | 'firstline'
%   opts.outer_frac    = 0.20    % fraction of edge width in kx/ky to use
%   opts.max_pages     = 32      % pages to sample for covariance estimate
%   opts.max_samples   = 20000   % cap number of samples used
%   opts.robust_clip   = true    % MAD clip per channel
%   opts.clip_k        = 6       % MAD multiplier
%   opts.lambda        = 1e-6    % regularisation
%   opts.seed          = 0       % RNG seed for reproducibility (0 = don't set)
%
% info:
%   .psi, .L, .p, .nSamples, .method

    if nargin < 3 || isempty(opts), opts = struct(); end
    opts = set_defaults(opts);

    origSize = size(kspace);
    if numel(origSize) < 3
        kspace_whitened = kspace;
        info = struct('psi',[],'L',[],'p',0,'nSamples',0,'method','none');
        return;
    end

    Nx = origSize(1);
    Ny = origSize(2);
    Nc = origSize(end);

    if Nc <= 1
        kspace_whitened = kspace;
        info = struct('psi',[],'L',[],'p',0,'nSamples',0,'method','singleChannel');
        return;
    end

    nTotal = numel(kspace);
    denom  = Nx * Ny * Nc;
    if mod(nTotal, denom) ~= 0
        error('noise_whiten:BadShape', ...
            'numel(kspace)=%d not divisible by Nx*Ny*Nc=%d. Channels must be last dim.', ...
            nTotal, denom);
    end
    Npages = nTotal / denom;

    % View as [Nx Ny Npages Nc]
    k4 = reshape(kspace, [Nx, Ny, Npages, Nc]);

    % ---------- Robust noise sample collection ----------
    eta = collect_noise_samples(k4, opts);   % [Nc x nSamples]

    % demean
    eta = eta - mean(eta, 2);

    % optional robust clipping (reduces spikes / residual signal influence)
    if opts.robust_clip
        eta = mad_clip(eta, opts.clip_k);
    end

    nSamples = size(eta, 2);
    if nSamples < 2
        kspace_whitened = kspace;
        info = struct('psi',[],'L',[],'p',0,'nSamples',nSamples,'method',opts.method);
        return;
    end

    % ---------- Covariance ----------
    psi = (eta * eta') / (nSamples - 1);

    % Regularise (scale by average variance)
    lam = opts.lambda * real(trace(psi))/Nc;
    if ~isfinite(lam) || lam <= 0, lam = opts.lambda; end
    psi = psi + lam * eye(Nc, 'like', psi);

    [L, p] = chol(psi, 'lower');

    % ---------- Apply whitening to all data ----------
    X = reshape(permute(k4, [4 1 2 3]), Nc, []); % [Nc x (Nx*Ny*Npages)]

    if p == 0
        Xw = L \ X;
    else
        % fallback: eigen whitening (handles non-PD)
        psiH = (psi + psi')/2;
        [V, D] = eig(psiH);
        d = real(diag(D));
        d(d <= 0) = eps(class(d));
        Winv = diag(1 ./ sqrt(d)) * V';
        Xw = Winv * X;
        L = [];
    end

    k4w = permute(reshape(Xw, [Nc, Nx, Ny, Npages]), [2 3 4 1]);
    kspace_whitened = reshape(k4w, origSize);

    info = struct('psi',psi,'L',L,'p',p,'nSamples',nSamples,'method',opts.method);
end


% ================= helper functions =================

function opts = set_defaults(opts)
    def.method      = 'outer';     % 'outer' | 'firstline'
    def.outer_frac  = 0.20;
    def.max_pages   = 32;
    def.max_samples = 20000;
    def.robust_clip = true;
    def.clip_k      = 6;
    def.lambda      = 1e-6;
    def.seed        = 0;

    f = fieldnames(def);
    for i = 1:numel(f)
        if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
            opts.(f{i}) = def.(f{i});
        end
    end
end

function eta = collect_noise_samples(k4, opts)
% k4: [Nx Ny Npages Nc]
    [Nx, Ny, Npages, Nc] = size(k4);

    if opts.seed ~= 0
        rng(opts.seed);
    end

    % Choose pages to sample (spread out; random is fine too)
    np = min(opts.max_pages, Npages);
    if np == Npages
        pages = 1:Npages;
    else
        pages = unique(round(linspace(1, Npages, np)));
    end

    switch lower(opts.method)
        case 'firstline'
            % Use kx=1, all ky, across selected pages: [Nc x (Ny*np)]
            tmp = k4(1, :, pages, :);              % [1 Ny np Nc]
            tmp = reshape(tmp, [], Nc);            % [(Ny*np) x Nc]
            eta = tmp.';                           % [Nc x nSamples]

        otherwise % 'outer'
            % Outer ring mask in kx/ky (noise-dominant in most cases)
            bx = max(1, floor(opts.outer_frac * Nx));
            by = max(1, floor(opts.outer_frac * Ny));

            mask = false(Nx, Ny);
            mask(1:bx, :) = true;
            mask(end-bx+1:end, :) = true;
            mask(:, 1:by) = true;
            mask(:, end-by+1:end) = true;

            % Collect samples from masked region across selected pages
            nMasked = nnz(mask);
            nPerPage = nMasked;
            nTot = nPerPage * numel(pages);

            % Reshape each selected page to [Nx*Ny x Nc] and take masked rows
            eta = zeros(Nc, min(nTot, opts.max_samples), 'like', k4);

            writePos = 1;
            maskv = mask(:);

            for ii = 1:numel(pages)
                p = pages(ii);
                page = reshape(k4(:,:,p,:), [], Nc);     % [Nx*Ny x Nc]
                samp = page(maskv, :);                   % [nMasked x Nc]

                % Optional random subsample to cap size
                if nTot > opts.max_samples
                    % take a fair share from each page
                    take = max(1, floor(opts.max_samples / numel(pages)));
                    if size(samp,1) > take
                        idx = randperm(size(samp,1), take);
                        samp = samp(idx, :);
                    end
                end

                nThis = size(samp, 1);
                if writePos + nThis - 1 > size(eta,2)
                    nThis = size(eta,2) - writePos + 1;
                    samp = samp(1:nThis, :);
                end

                eta(:, writePos:writePos+nThis-1) = samp.';  % [Nc x nThis]
                writePos = writePos + nThis;

                if writePos > size(eta,2)
                    break;
                end
            end

            eta = eta(:, 1:writePos-1);
    end
end

function Xc = mad_clip(X, k)
% Clip magnitude outliers per channel using MAD on |X|
% X: [Nc x nSamples]
    mag = abs(X);
    med = median(mag, 2);
    madv = median(abs(mag - med), 2) + eps(class(mag));
    thr = med + k * 1.4826 * madv;   % ~sigma for Gaussian

    % scale down samples exceeding threshold (preserve phase)
    scale = min(1, thr ./ (mag + eps(class(mag))));
    Xc = X .* scale;
end
