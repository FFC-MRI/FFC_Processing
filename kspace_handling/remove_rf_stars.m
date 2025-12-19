function [k_clean, dbg] = remove_rf_stars(k, opts)
%SUBTRACT_READOUT_TONES_LS_SAFE Robust multi-tone subtraction along readout without blanking.
%
% k dims: [Nx Ny Nz Nt ... Nc] (readout, phase, slices, times, ..., coils)
%
% Key improvements vs naive LS:
%  - robust weights downweight real MR signal in the fit region
%  - only apply subtraction if model matches background sufficiently (gate)
%  - cap subtraction energy relative to background energy (prevents over-subtract)
%
% opts:
%   .nTones            = 3
%   .outerFrac         = 0.35
%   .excludeCenterFrac = 0.35
%   .excludeDCBins     = 3
%   .minPromDB         = 8
%   .lambda            = 1e-6
%   .robustC           = 3.0      % smaller => stronger downweighting of big samples
%   .minCorr           = 0.5      % gate: require correlation(model, y) on fit region
%   .maxSubFrac        = 0.8      % cap: ||model|| <= maxSubFrac*||y|| (on fit region)
%   .verbose           = false

    if nargin < 2, opts = struct(); end
    opts = fill_defaults(opts);

    sz = size(k);
    Nx = sz(1);

    % coils last if present
    if numel(sz) >= 5
        Nc = sz(end);
        hasCoils = true;
    else
        Nc = 1;
        hasCoils = false;
    end

    % reshape to [Nx nLines Nc]
    if hasCoils
        k3 = reshape(k, Nx, [], Nc);
    else
        k3 = reshape(k, Nx, [], 1);
    end
    nLines = size(k3,2);

    % 1D fit mask along readout
    fitMask = bg_mask_1d(Nx, opts.outerFrac, opts.excludeCenterFrac);
    nb = nnz(fitMask);
    if nb < 24
        error('Fit mask too small. Increase outerFrac or reduce excludeCenterFrac.');
    end

    n = (0:Nx-1).';

    % ---------- Global tone frequency detection on pooled robust background ----------
    pooled = pooled_background(k3, fitMask);
    [wList, specDB, fgrid] = pick_tones_from_pooled(pooled, opts);

    if isempty(wList)
        k_clean = k;
        dbg = struct('w',[],'fitMask',fitMask,'specDB',specDB,'fgrid',fgrid);
        return;
    end

    K = numel(wList);

    % Precompute basis matrices
    E  = exp(1j * (n * wList(:).'));     % [Nx x K]
    Eb = E(fitMask, :);                 % [nb x K]

    k3_clean = k3;
    A = zeros(K, nLines, Nc, 'like', k3);
    applied = false(nLines, Nc);

    for c = 1:Nc
        for L = 1:nLines
            y = k3(:,L,c);
            yb = y(fitMask);

            % ---- robust weights: downweight large |yb| (likely object signal leakage) ----
            % scale from MAD of |yb|
            s = median(abs(abs(yb) - median(abs(yb)))) * 1.4826;
            if ~isfinite(s) || s <= 0, s = max(abs(yb)) + eps; end
            r = abs(yb) / (opts.robustC*s + eps);
            w = 1 ./ (1 + r.^2);                     % Cauchy-like weights in (0,1]
            W = w(:);

            % ---- weighted ridge LS solve: a = (Eb'*W*Eb + lam I)\(Eb'*W*yb) ----
            % implement via scaling rows
            Ews = Eb .* W;                            % [nb x K]
            R = (Eb' * Ews) + opts.lambda * eye(K, 'like', Eb);
            b = (Eb' * (W .* yb));
            a = R \ b;

            yhat_b = Eb * a;

            % ---- gate: only subtract if it really matches background ----
            corrVal = abs( (yhat_b' * yb) / (norm(yhat_b)*norm(yb) + eps) );
            if corrVal < opts.minCorr
                continue;
            end

            % ---- cap: don't subtract more energy than a fraction of background energy ----
            subFrac = norm(yhat_b) / (norm(yb) + eps);
            if subFrac > opts.maxSubFrac
                a = a * (opts.maxSubFrac / subFrac);
                yhat_b = Eb * a;
            end

            % apply subtraction to full line
            yhat = E * a;
            k3_clean(:,L,c) = y - yhat;

            A(:,L,c) = a;
            applied(L,c) = true;
        end
    end

    k_clean = reshape(k3_clean, sz);

    dbg = struct();
    dbg.w = wList;
    dbg.fitMask = fitMask;
    dbg.specDB = specDB;
    dbg.fgrid = fgrid;
    dbg.A = A;
    dbg.appliedFrac = nnz(applied) / numel(applied);
end

% ---------------- helpers ----------------

function opts = fill_defaults(opts)
    def.nTones            = 3;
    def.outerFrac         = 0.35;
    def.excludeCenterFrac = 0.35;
    def.excludeDCBins     = 3;
    def.minPromDB         = 8;
    def.lambda            = 1e-6;
    def.robustC           = 3.0;
    def.minCorr           = 0.5;
    def.maxSubFrac        = 0.8;
    def.verbose           = false;

    f = fieldnames(def);
    for i=1:numel(f)
        if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
    end
end

function m = bg_mask_1d(Nx, outerFrac, excludeCenterFrac)
    m = false(Nx,1);
    nOuter = max(1, round(outerFrac*Nx));
    m(1:nOuter) = true;
    m(end-nOuter+1:end) = true;

    c = floor(Nx/2)+1;
    hw = round(excludeCenterFrac*Nx/2);
    m(max(1,c-hw):min(Nx,c+hw)) = false;
end

function pooled = pooled_background(k3, fitMask)
    % k3: [Nx nLines Nc]
    kb = k3(fitMask,:,:);               % [nb nLines Nc]
    kb = reshape(kb, size(kb,1), []);   % [nb (nLines*Nc)]
    pooled = median(kb, 2);
    pooled = pooled - mean(pooled);
end

function [wList, specDB, fgrid] = pick_tones_from_pooled(pooled, opts)
    nb = numel(pooled);
    nfft = 2^nextpow2(max(2048, 8*nb));
    X = fftshift(fft(pooled, nfft));
    mag = abs(X);
    specDB = 20*log10(mag + eps);
    specDB = specDB - max(specDB);
    fgrid = linspace(-pi, pi, nfft).';

    mid = floor(nfft/2)+1;
    specDB(mid-opts.excludeDCBins:mid+opts.excludeDCBins) = -Inf;

    % pick peaks and always include +/- pairs
    work = specDB;
    finiteVals = work(isfinite(work));
    if isempty(finiteVals)
        wList = [];
        return;
    end
    floorDB = median(finiteVals);

    wList = [];
    for t = 1:opts.nTones
        [pk, idx] = max(work);
        if ~isfinite(pk), break; end
        prom = pk - floorDB;
        if prom < opts.minPromDB, break; end

        w = fgrid(idx);
        wList(end+1,1) = w; %#ok<AGROW>

        % add conjugate partner
        [~, idx2] = min(abs(fgrid + w));
        if idx2 ~= idx
            wList(end+1,1) = fgrid(idx2); %#ok<AGROW>
        end

        sup = max(3, round(0.01*numel(fgrid)));
        work(max(1,idx-sup):min(numel(fgrid),idx+sup)) = -Inf;
        work(max(1,idx2-sup):min(numel(fgrid),idx2+sup)) = -Inf;
    end
end
