function [k_clean, model, dbg] = remove_rf_tones_kspace(k, opts)
%REMOVE_RF_TONES_KSPACE Fit & subtract multiple complex sinusoidal interferers.
%
%   k_clean = remove_rf_tones_kspace(k, opts)
%
% k: complex k-space, at minimum [Nx Ny]. Can be higher-D; will be reshaped
%    so that dim1 = readout samples (kx) and dim2 = "lines" (everything else).
%
% This estimates tone frequencies from "background" samples (outer kx),
% then fits complex exponentials exp(1j*w*n) and subtracts them.
%
% opts fields (all optional):
%   .nTones           = number of tones to fit (default 3)
%   .outerFrac        = fraction of outer kx used as background (default 0.25)
%   .excludeCenterFrac= fraction around kx center excluded from background (default 0.20)
%   .perLineAmpPhase  = fit separate complex amplitudes per line (default true)
%   .perLineFreqRefine= small per-line freq refinement (default false)
%   .refineBins       = +/- bins for refinement (default 2)
%   .minPeakProm      = min peak prominence in spectrum (default 6) in dB
%   .maxIter          = max iterations for tone selection (default nTones)
%   .detrendOrder     = remove polynomial baseline on bg samples (default 0)
%   .lambda           = ridge regularization for LS (default 1e-8)
%   .verbose          = true/false (default false)
%
% Outputs:
%   k_clean: cleaned k-space, same size as k
%   model  : fitted interference model in reshaped form [Nx nLines]
%   dbg    : diagnostics (freqs, spectrum, mask, etc.)

    if nargin < 2, opts = struct(); end
    opts = fill_defaults(opts);

    sz = size(k);
    Nx = sz(1);
    k2 = reshape(k, Nx, []);     % [Nx x nLines]
    nLines = size(k2, 2);

    n = (0:Nx-1).';              % sample index (readout time up to scaling)

    % --- Build "background" mask over kx (outer samples, exclude center) ---
    bgMask = background_mask(Nx, opts.outerFrac, opts.excludeCenterFrac);
    nb = nnz(bgMask);
    if nb < max(16, 0.05*Nx)
        error('Background mask too small; adjust outerFrac/excludeCenterFrac.');
    end

    % --- Build a pooled background signal for frequency detection ---
    % Robust average over lines: median across lines to suppress object signal
    kb = k2(bgMask, :);                       % [nb x nLines]
    pooled = median(kb, 2);                   % [nb x 1]
    pooled = pooled - mean(pooled);           % remove DC

    % Optional polynomial detrend on pooled
    if opts.detrendOrder > 0
        x = linspace(-1, 1, nb).';
        p = polyfit(x, pooled, opts.detrendOrder);
        pooled = pooled - polyval(p, x);
    end

    % --- Estimate tone frequencies iteratively from pooled spectrum ---
    [toneBins, toneFreqsRad, specDB, fgrid] = pick_tones(pooled, opts);

    if opts.verbose
        fprintf('Selected %d tone(s).\n', numel(toneBins));
        disp(toneFreqsRad(:).');
    end

    % --- Fit & subtract tones ---
    % Design matrix for tones at the chosen frequencies
    W = toneFreqsRad(:).';                 % [1 x K] rad/sample
    K = numel(W);
    E = exp(1j * (n * W));                 % [Nx x K]

    model2 = zeros(Nx, nLines, 'like', k2);

    if opts.perLineAmpPhase
        % Fit amplitudes per line using background samples only
        Eb = E(bgMask, :);                 % [nb x K]
        % Regularized normal equations: a = (Eb'*Eb + lam I)\(Eb'*y)
        R = Eb' * Eb + opts.lambda * eye(K, 'like', Eb);
        Rt = (R \ Eb');                    % [K x nb] (precompute)
        for L = 1:nLines
            y = k2(bgMask, L);
            a = Rt * y;                    % [K x 1]
            model2(:, L) = E * a;
        end

    else
        % One global amplitude vector across all lines (less flexible)
        Eb = E(bgMask, :);
        y  = pooled;
        a  = (Eb' * Eb + opts.lambda * eye(K, 'like', Eb)) \ (Eb' * y);
        model2 = E * a;
        model2 = repmat(model2, 1, nLines);
    end

    % Optional: small per-line frequency refinement around selected bins
    if opts.perLineFreqRefine && K > 0
        model2 = refine_per_line(k2, model2, bgMask, n, toneBins, opts);
    end

    k2_clean = k2 - model2;

    % reshape back
    k_clean = reshape(k2_clean, sz);
    model   = reshape(model2, sz);

    dbg = struct();
    dbg.bgMask = bgMask;
    dbg.toneBins = toneBins;
    dbg.toneFreqsRadPerSample = toneFreqsRad;
    dbg.spectrumDB = specDB;
    dbg.fgrid = fgrid;
end

% ---------------- helpers ----------------

function opts = fill_defaults(opts)
    def.nTones            = 3;
    def.outerFrac         = 0.25;
    def.excludeCenterFrac = 0.20;
    def.perLineAmpPhase   = true;
    def.perLineFreqRefine = false;
    def.refineBins        = 2;
    def.minPeakProm       = 6;     % dB
    def.maxIter           = [];    % default = nTones
    def.detrendOrder      = 0;
    def.lambda            = 1e-8;
    def.verbose           = false;

    f = fieldnames(def);
    for i=1:numel(f)
        if ~isfield(opts, f{i}) || isempty(opts.(f{i}))
            opts.(f{i}) = def.(f{i});
        end
    end
    if isempty(opts.maxIter), opts.maxIter = opts.nTones; end
end

function m = background_mask(Nx, outerFrac, excludeCenterFrac)
    m = false(Nx,1);
    nOuter = max(1, round(outerFrac * Nx));
    m(1:nOuter) = true;
    m(end-nOuter+1:end) = true;

    % exclude a window around center (often contains object signal leakage)
    c = floor(Nx/2)+1;
    hw = round(excludeCenterFrac * Nx / 2);
    lo = max(1, c-hw);
    hi = min(Nx, c+hw);
    m(lo:hi) = false;
end

function [toneBins, toneFreqsRad, specDB, fgrid] = pick_tones(pooled, opts)
    % FFT on pooled background vector. Use zero-padding for finer binning.
    nb = numel(pooled);
    nfft = 2^nextpow2(max(1024, 4*nb));
    X = fftshift(fft(pooled, nfft));
    mag = abs(X);

    % Avoid selecting DC
    mid = floor(nfft/2)+1;
    mag(mid-1:mid+1) = 0;

    % Convert to dB for peak picking
    specDB = 20*log10(mag + eps);
    specDB = specDB - max(specDB);  % normalize (0 dB max)
    fgrid = linspace(-pi, pi, nfft).';  % rad/sample, aligned with fftshift

    toneBins = [];
    toneFreqsRad = [];

    resid = mag;
    for it = 1:opts.maxIter
        [pk, idx] = max(resid);
        if isempty(pk) || pk <= 0
            break;
        end

        % Check peak prominence-ish: compare to median floor
        floorVal = median(resid(resid > 0));
        promDB = 20*log10((pk+eps)/(floorVal+eps));
        if promDB < opts.minPeakProm
            break;
        end

        toneBins(end+1) = idx; %#ok<AGROW>
        toneFreqsRad(end+1) = fgrid(idx); %#ok<AGROW>

        % Notch out a small neighborhood around the selected peak so we can find others
        notch = max(3, round(0.01*nfft));
        lo = max(1, idx-notch);
        hi = min(nfft, idx+notch);
        resid(lo:hi) = 0;
    end
end

function model2 = refine_per_line(k2, model2, bgMask, n, toneBins, opts)
    % Optional refinement: allow each line to shift each tone by a few FFT bins.
    % This can help if the interferer drifts slightly over time/lines.
    % Note: this is heavier; keep refineBins small.

    [Nx, nLines] = size(k2);
    K = numel(toneBins);
    if K == 0, return; end

    % For refinement we build candidate frequencies around each bin
    % by converting bin index offsets to rad/sample on an implicit grid.
    % We need an approximate nfft; infer from bin indices spacing:
    % We'll approximate using nfft = max bin index range * 2.
    nfft = max(toneBins) * 2; % crude, but OK for small +/- bin refinement
    fgrid = linspace(-pi, pi, nfft).'; 

    for L = 1:nLines
        yb = k2(bgMask, L);
        bestModel = zeros(Nx,1,'like',k2);
        for k = 1:K
            idx0 = toneBins(k);
            candIdx = (idx0-opts.refineBins):(idx0+opts.refineBins);
            candIdx = candIdx(candIdx>=1 & candIdx<=nfft);
            bestErr = inf;
            bestMk = zeros(Nx,1,'like',k2);

            for idx = candIdx
                w = fgrid(idx);
                e = exp(1j*w*n);
                eb = e(bgMask);

                a = (eb' * eb + opts.lambda) \ (eb' * yb);
                mk = e * a;
                err = norm((yb - mk(bgMask)), 2);

                if err < bestErr
                    bestErr = err;
                    bestMk = mk;
                end
            end
            bestModel = bestModel + bestMk;
            % Update residual for next tone (greedy)
            yb = yb - bestMk(bgMask);
        end
        model2(:,L) = bestModel;
    end
end
