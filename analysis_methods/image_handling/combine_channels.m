function [combined_images] = combine_channels(images, opts, obj)
%COMBINE_CHANNELS  Multicoil reconstruction with per-frame masked self-calibrated adaptive combine
%
% opts is a logical coil mask (dimension 6 is coils).
% If a single coil is selected, returns that coil (debug mode).
%
% Assumes noise-whitening has already been applied upstream.

    % Select requested coils (dim 6 is coil dimension)
    images = images(:,:,:,:,:,opts);

    tempdim = size(images);

    % Reduce everything to 4 dimensions: (X,Y,frames,nCoils)
    images_temp = reshape(images, tempdim(1), tempdim(2), [], size(images,6));

    if length(opts) > 1
        n_channels = size(images_temp, 4);

        if n_channels > 1

            % ----------------------------
            % Settings (optionally from obj)
            % ----------------------------
            sigma = 7; % good default for 128x128 in vivo
            if isstruct(obj) && isfield(obj, 'coil_sens_sigma') && ~isempty(obj.coil_sens_sigma)
                sigma = obj.coil_sens_sigma;
            end

            maskFrac = 0.06; % fraction of 95th percentile of RSS
            if isstruct(obj) && isfield(obj, 'coil_sens_mask_frac') && ~isempty(obj.coil_sens_mask_frac)
                maskFrac = obj.coil_sens_mask_frac;
            end

            lambda = 0; % optional extra regularisation in denominator
            if isstruct(obj) && isfield(obj, 'coil_sens_lambda') && ~isempty(obj.coil_sens_lambda)
                lambda = obj.coil_sens_lambda;
            end

            nFrames = size(images_temp, 3);
            combined = zeros(tempdim(1), tempdim(2), nFrames, 'like', images_temp);

            % ----------------------------
            % Per-frame sensitivity estimation + combine (no global ref frame)
            % ----------------------------
            for k = 1:nFrames
                Ik = images_temp(:,:,k,:);                   % (X,Y,1,coil)
                rss = sqrt(sum(abs(Ik).^2, 4));
                rss = rss + eps(class(rss));

                % Robust object mask from RSS
                p95 = prctile(rss(:), 95);
                thr = max(eps(class(rss)), maskFrac * p95);
                objmask = rss > thr;

                % Smooth mask (soft weighting), then do normalised convolution smoothing
                objmask_f = imgaussfilt(single(objmask), max(1, sigma/2), ...
                    'FilterSize', max(3, 2*ceil(3*max(1,sigma/2))+1));

                coil_norm = Ik ./ rss;                       % (X,Y,1,coil)
                sens = local_masked_smooth_complex_2d(coil_norm, objmask_f, sigma); % (X,Y,1,coil)

                % Normalise sens so sum(|sens|^2) ~ 1
                sens_norm = sqrt(sum(abs(sens).^2, 4));
                sens = sens ./ (sens_norm + eps(class(sens_norm)));

                % Roemer / adaptive complex combination
                num = sum(conj(sens) .* Ik, 4);
                den = sum(abs(sens).^2, 4) + lambda + eps(class(num));
                combined(:,:,k) = num ./ den;                % complex
            end

            % Reshape back to expected dimensionality
            combined_images = reshape(combined, ...
                [tempdim(1), tempdim(2), obj.slices, obj.n_timepoints, obj.n_fieldpoints]);

        else
            % Only one channel selected
            combined_images = reshape(images_temp(:,:,:,1), ...
                [tempdim(1), tempdim(2), obj.slices, obj.n_timepoints, obj.n_fieldpoints]);
        end

    else
        % Displays the individual channel data for debugging purposes
        combined_images = images(:,:,:,:,:,1);
    end
end

% ---- helper: masked, normalised smoothing of complex 2D fields (per frame) ----
function out = local_masked_smooth_complex_2d(in, mask_f, sigma)
    % in: (X,Y,1,nCoils) complex
    % mask_f: (X,Y) soft mask in [0..1]
    % out: (X,Y,1,nCoils) complex

    out = zeros(size(in), 'like', in);

    denom = imgaussfilt(mask_f, sigma, 'FilterSize', max(3, 2*ceil(3*sigma)+1));
    denom = denom + eps(class(denom));

    nCoils = size(in,4);
    for c = 1:nCoils
        re_num = imgaussfilt(real(in(:,:,1,c)) .* mask_f, sigma, 'FilterSize', max(3, 2*ceil(3*sigma)+1));
        im_num = imgaussfilt(imag(in(:,:,1,c)) .* mask_f, sigma, 'FilterSize', max(3, 2*ceil(3*sigma)+1));
        out(:,:,1,c) = complex(re_num ./ denom, im_num ./ denom);
    end
end
