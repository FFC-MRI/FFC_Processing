function combined = combine_adaptive_smoothed(images, opts)
%COMBINE_ADAPTIVE_SMOOTHED_REFAWARE
%
% Adaptive complex coil combine for whitened physical/virtual channels.
%
% images: [X Y slices time field coils] complex
% opts: either:
%   - numeric vector of selected coil indices, e.g. 1:nCoils
%   - struct with fields:
%       opts.coils              = selected coil indices
%       opts.sigma              = smoothing sigma, default 7
%       opts.maskFrac           = object mask fraction, default 0.06
%       opts.referenceMode      = 'rss' | 'signalmask' | 'signalcoil'
%
% referenceMode:
%
%   'rss'
%       Original behaviour. Use RSS for both mask and normalisation.
%
%   'signalmask'
%       Use final channel, assumed to be the true NMR coil, to define the
%       object mask. Use RSS for coil normalisation.
%
%   'signalcoil'
%       Use final channel for both mask and normalisation.
%
% Output:
%   combined: [X Y slices time field] complex

    % ------------------------------------------------------------
    % Backwards-compatible opts handling
    % ------------------------------------------------------------

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    if isnumeric(opts) || islogical(opts)
        tmp = struct();
        tmp.coils = opts;
        opts = tmp;
    end

    if ~isfield(opts, 'coils') || isempty(opts.coils)
        opts.coils = 1:size(images, 6);
    end

    if ~isfield(opts, 'sigma') || isempty(opts.sigma)
        opts.sigma = 7;
    end

    if ~isfield(opts, 'maskFrac') || isempty(opts.maskFrac)
        opts.maskFrac = 0.06;
    end

    if ~isfield(opts, 'referenceMode') || isempty(opts.referenceMode)
        opts.referenceMode = 'rss';
    end

    coils = opts.coils;
    sigma = opts.sigma;
    maskFrac = opts.maskFrac;
    referenceMode = lower(char(opts.referenceMode));

    % ------------------------------------------------------------
    % Select coils
    % ------------------------------------------------------------

    images = images(:,:,:,:,:,coils);

    sz = size(images);
    sz(end+1:6) = 1;

    X = sz(1);
    Y = sz(2);
    S = sz(3);
    T = sz(4);
    F = sz(5);
    C = sz(6);

    tmp = reshape(images, X, Y, S*T*F, C);
    combined_tmp = zeros(X, Y, S*T*F, 'like', images);

    % ------------------------------------------------------------
    % Adaptive combine
    % ------------------------------------------------------------

    for k = 1:size(tmp,3)

        Ik = tmp(:,:,k,:);  % [X Y 1 C]

        rss = sqrt(sum(abs(Ik).^2, 4)) + eps('single');
        sig = abs(Ik(:,:,1,end)) + eps('single');

        switch referenceMode

            case 'rss'
                % Original behaviour
                maskBase = rss;
                normBase = rss;

            case 'signalmask'
                % Use true signal coil for mask, but RSS for normalisation
                maskBase = sig;
                normBase = rss;

            case 'signalcoil'
                % Use true signal coil for both mask and normalisation
                maskBase = sig;
                normBase = sig;

            otherwise
                error('Unknown opts.referenceMode: %s. Use ''rss'', ''signalmask'', or ''signalcoil''.', ...
                    referenceMode);
        end

        p95 = prctile(maskBase(:), 95);
        thr = max(eps('single'), maskFrac * p95);
        mask = maskBase > thr;

        maskSigma = max(1, sigma/2);
        maskFiltSz = max(3, 2*ceil(3*maskSigma)+1);

        mask_f = imgaussfilt(single(mask), maskSigma, ...
            'FilterSize', maskFiltSz);

        coil_norm = Ik ./ normBase;

        sens = local_masked_smooth_complex_2d(coil_norm, mask_f, sigma);

        sens_norm = sqrt(sum(abs(sens).^2, 4));
        sens = sens ./ (sens_norm + eps('single'));

        num = sum(conj(sens) .* Ik, 4);
        den = sum(abs(sens).^2, 4) + eps('single');

        combined_tmp(:,:,k) = num ./ den;
    end

    combined = reshape(combined_tmp, X, Y, S, T, F);
end


function out = local_masked_smooth_complex_2d(in, mask_f, sigma)
%LOCAL_MASKED_SMOOTH_COMPLEX_2D
%
% in: [X Y 1 C]

    out = zeros(size(in), 'like', in);

    filtSz = max(3, 2*ceil(3*sigma)+1);
    denom = imgaussfilt(mask_f, sigma, 'FilterSize', filtSz) + eps('single');

    C = size(in, 4);

    for c = 1:C
        re_num = imgaussfilt(real(in(:,:,1,c)) .* mask_f, sigma, ...
            'FilterSize', filtSz);

        im_num = imgaussfilt(imag(in(:,:,1,c)) .* mask_f, sigma, ...
            'FilterSize', filtSz);

        out(:,:,1,c) = complex(re_num ./ denom, im_num ./ denom);
    end
end