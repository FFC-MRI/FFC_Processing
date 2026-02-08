function mag_corr = bias_correct_n4lite(mag, obj)
%BIAS_CORRECT_N4LITE  Smooth multiplicative bias-field correction on magnitude images
%
% mag:     [X Y slices time field] magnitude
% output:  same size, bias-corrected magnitude
%
% Robust to motion/phase drift because it operates on magnitude only.

    sz = size(mag);
    sz(end+1:5) = 1;
    X = sz(1); Y = sz(2); S = sz(3); T = sz(4); F = sz(5);

    % ----- Tunables (good starting points for 128x128 in vivo) -----
    sigma_field = 10;         % pixels; smoothness of bias field (try 16–28 for 128x128)
    maskFrac    = 0.14;       % threshold fraction of p95 within each frame
    epsFloor    = 1e-8;       % to avoid log(0)

    if isfield(obj,'bias_sigma') && ~isempty(obj.bias_sigma), sigma_field = obj.bias_sigma; end
    if isfield(obj,'bias_mask_frac') && ~isempty(obj.bias_mask_frac), maskFrac = obj.bias_mask_frac; end
    if isfield(obj,'bias_eps') && ~isempty(obj.bias_eps), epsFloor = obj.bias_eps; end

    filtSz = max(3, 2*ceil(3*sigma_field)+1);

    mag_corr = zeros(X,Y,S,T,F, 'like', mag);

    for s = 1:S
        for t = 1:T
            for f = 1:F
                I = mag(:,:,s,t,f);

                % --- Build a robust soft mask from this frame ---
                p95 = prctile(I(:), 95);
                thr = max(eps(class(I)), maskFrac * p95);
                mask = I > thr;

                % soften the mask so the field doesn't ring at edges
                mask_f = imgaussfilt(single(mask), max(1, sigma_field/3), ...
                    'FilterSize', max(3, 2*ceil(3*max(1,sigma_field/3))+1));

                % --- Estimate multiplicative bias field in log domain ---
                % logI = log( true * bias ) = log(true) + log(bias)
                logI = log(double(I) + epsFloor);

                % Masked, normalised smoothing (prevents background dominating)
                num = imgaussfilt(logI .* mask_f, sigma_field, 'FilterSize', filtSz);
                den = imgaussfilt(mask_f,        sigma_field, 'FilterSize', filtSz) + eps('single');
                logBias = num ./ den;

                bias = exp(logBias);  % smooth multiplicative field

                % --- Correct and rescale (optional) ---
                Icorr = double(I) ./ (bias + epsFloor);

                % Rescale so median in mask is preserved (prevents arbitrary gain changes)
                med_before = median(double(I(mask)));
                med_after  = median(Icorr(mask));
                if isfinite(med_before) && isfinite(med_after) && med_after > 0
                    Icorr = Icorr * (med_before / med_after);
                end

                mag_corr(:,:,s,t,f) = cast(Icorr, 'like', mag);
            end
        end
    end
end
