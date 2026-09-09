function [mag_corr, bias, mask3] = map_guided_bias_correct(mag, maps, obj)
%MAP_GUIDED_BIAS_CORRECT  Coil-sensitivity-guided bias correction (with optional upgrade).
%
% Upgrade:
%   Uses a normalised sensitivity RSS field for bias estimation:
%       norm = sqrt(sum_c |S|^2)
%       srN  = norm / median(norm(mask))
%       bias = smooth(srN)   (then re-normalise to median=1)
%
% This makes the bias field less sensitive to residual scaling differences
% and reduces any anatomy imprinting indirectly introduced via map estimation.

    % ---- shapes ----
    sz = size(mag);  sz(end+1:5) = 1;
    X = sz(1); Y = sz(2); S = sz(3); T = sz(4); F = sz(5);

    mapsz = size(maps); mapsz(end+1:6) = 1;
    if mapsz(1) ~= X || mapsz(2) ~= Y || mapsz(3) ~= S
        error('map_guided_bias_correct:sizeMismatch', ...
            'mag is [%d %d %d ...], maps is [%d %d %d ...].', X,Y,S, mapsz(1),mapsz(2),mapsz(3));
    end

    % ---- tunables ----
    sigma    = getprop_default(obj, 'map_bias_sigma', 18);
    maskFrac = getprop_default(obj, 'map_bias_mask_frac', 0.06);
    bmin     = getprop_default(obj, 'map_bias_clip_min', 0.4);
    bmax     = getprop_default(obj, 'map_bias_clip_max', 2.5);
    dilr     = getprop_default(obj, 'map_bias_dilate', 2);

    filtSz = max(3, 2*ceil(3*sigma)+1);

    % ---- sensitivity norm field (coil-driven) ----
    % norm = sqrt(sum |S|^2)  (this is the upgrade base)
    norm = sqrt(sum(abs(maps).^2, 6));         % [X Y S 1 1]
    norm = reshape(norm, [X Y S]);             % [X Y S]
    norm = single(norm) + eps('single');

    bias  = ones(X,Y,S,'single');
    mask3 = false(X,Y,S);

    for s = 1:S
        sr = norm(:,:,s);

        % Mask from norm (coil-driven)
        p95 = prctile(sr(:), 95);
        thr = max(eps('single'), maskFrac * p95);
        m = sr > thr;

        % Clean mask (fill + dilate)
        m = imfill(m, 'holes');
        if dilr > 0
            m = imdilate(m, strel('disk', dilr));
        end

        % ---- OPTIONAL UPGRADE STEP ----
        % Normalise sr inside the mask so it has median 1 before smoothing.
        % This removes arbitrary scaling and focuses smoothing on spatial variation.
        med_sr = median(sr(m));
        if ~isfinite(med_sr) || med_sr <= 0, med_sr = 1; end
        srN = sr / med_sr;

        % Smooth the normalised sensitivity norm to get low-frequency bias
        b = imgaussfilt(srN, sigma, 'FilterSize', filtSz);

        % Re-normalise bias to median 1 inside mask (safe)
        med_b = median(b(m));
        if ~isfinite(med_b) || med_b <= 0, med_b = 1; end
        b = b / med_b;

        % Clamp to avoid extreme amplification
        b = max(bmin, min(bmax, b));

        bias(:,:,s)  = b;
        mask3(:,:,s) = m;
    end

    % ---- apply only inside mask, leave background unchanged ----
    mag_corr = mag;
    for s = 1:S
        m = mask3(:,:,s);
        if ~any(m(:)), continue; end

        b = bias(:,:,s);
        for t = 1:T
            for f = 1:F
                I = mag(:,:,s,t,f);
                I(m) = I(m) ./ (b(m) + eps('single'));
                mag_corr(:,:,s,t,f) = I;
            end
        end
    end
end

function v = getprop_default(obj, name, default)
%GETPROP_DEFAULT Read class property if present, else default.
    if isprop(obj, name)
        v = obj.(name);
        if isempty(v)
            v = default;
        end
    else
        v = default;
    end
end
