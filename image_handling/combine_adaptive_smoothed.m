function combined = combine_adaptive_smoothed(images, opts, obj)
%COMBINE_ADAPTIVE_SMOOTHED  Self-calibrated adaptive combine via smoothed coil-normalised fields
%
% images: [X Y slices time field coils] complex (already trimmed to selected coils)
% opts:   1:nCoils (local index space)
% output: [X Y slices time field] complex

    images = images(:,:,:,:,:,opts);

    sigma = 7;
    maskFrac = 0.06;
    if isfield(obj,'coil_sens_sigma') && ~isempty(obj.coil_sens_sigma), sigma = obj.coil_sens_sigma; end
    if isfield(obj,'coil_sens_mask_frac') && ~isempty(obj.coil_sens_mask_frac), maskFrac = obj.coil_sens_mask_frac; end

    sz = size(images); sz(end+1:6) = 1;
    X=sz(1); Y=sz(2); S=sz(3); T=sz(4); F=sz(5); C=sz(6);

    tmp = reshape(images, X, Y, S*T*F, C);
    combined_tmp = zeros(X, Y, S*T*F, 'like', images);

    for k = 1:size(tmp,3)
        Ik = tmp(:,:,k,:);                           % (X,Y,1,C)
        rss = sqrt(sum(abs(Ik).^2,4)) + eps('single');

        p95 = prctile(rss(:),95);
        thr = max(eps('single'), maskFrac * p95);
        mask = rss > thr;

        mask_f = imgaussfilt(single(mask), max(1,sigma/2), ...
            'FilterSize', max(3, 2*ceil(3*max(1,sigma/2))+1));

        coil_norm = Ik ./ rss;
        sens = local_masked_smooth_complex_2d(coil_norm, mask_f, sigma);

        sens_norm = sqrt(sum(abs(sens).^2,4));
        sens = sens ./ (sens_norm + eps('single'));

        num = sum(conj(sens).*Ik,4);
        den = sum(abs(sens).^2,4) + eps('single');
        combined_tmp(:,:,k) = num ./ den;
    end

    combined = reshape(combined_tmp, X, Y, S, T, F);
end

function out = local_masked_smooth_complex_2d(in, mask_f, sigma)
% in: (X,Y,1,C)
    out = zeros(size(in), 'like', in);
    filtSz = max(3, 2*ceil(3*sigma)+1);
    denom = imgaussfilt(mask_f, sigma, 'FilterSize', filtSz) + eps('single');

    C = size(in,4);
    for c = 1:C
        re_num = imgaussfilt(real(in(:,:,1,c)).*mask_f, sigma, 'FilterSize', filtSz);
        im_num = imgaussfilt(imag(in(:,:,1,c)).*mask_f, sigma, 'FilterSize', filtSz);
        out(:,:,1,c) = complex(re_num./denom, im_num./denom);
    end
end
