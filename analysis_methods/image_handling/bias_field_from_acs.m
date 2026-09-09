function [bias, mask3] = bias_field_from_acs(kspace, opts, obj)
% Returns:
%   bias  [X Y S]  (median=1 inside mask)
%   mask3 [X Y S]  logical body mask (for safe application)

    kspace = kspace(:,:,:,:,:,opts);

    sz = size(kspace); sz(end+1:6)=1;
    X=sz(1); Y=sz(2); S=sz(3); T=sz(4); F=sz(5); C=sz(6);
    kspace = reshape(kspace, [X Y S T F C]);

    acs = 24;
    smoothSigma = 18;
    maskFrac = 0.45;
    if isfield(obj,'espirit_acs') && ~isempty(obj.espirit_acs), acs = obj.espirit_acs; end
    if isfield(obj,'bias_sigma') && ~isempty(obj.bias_sigma), smoothSigma = obj.bias_sigma; end
    if isfield(obj,'bias_mask_frac') && ~isempty(obj.bias_mask_frac), maskFrac = obj.bias_mask_frac; end

    cx = floor(X/2)+1; cy = floor(Y/2)+1;
    hx = floor(acs/2); hy = floor(acs/2);
    xr = max(1,cx-hx):min(X,cx+hx-1);
    yr = max(1,cy-hy):min(Y,cy+hy-1);

    filtSz = max(3, 2*ceil(3*smoothSigma)+1);

    bias  = ones(X,Y,S,'single');
    mask3 = false(X,Y,S);

    for s = 1:S
        % Pick best (t,f) by ACS energy (robust to drift/motion)
        E = zeros(T,F,'double');
        for t = 1:T
            for f = 1:F
                blk = kspace(xr,yr,s,t,f,:);      % [acs acs 1 1 1 C]
                E(t,f) = sum(abs(blk(:)).^2);
            end
        end
        [~,idx] = max(E(:));
        [tBest,fBest] = ind2sub([T F], idx);

        % ACS-only k-space for that best frame
        Kacs = zeros(X,Y,1,1,C,'like',kspace);
        Kacs(xr,yr,1,1,:) = reshape(kspace(xr,yr,s,tBest,fBest,:), [numel(xr) numel(yr) 1 1 C]);

        L = ifft2c(Kacs);                 % [X Y 1 1 C]
        L = reshape(L, [X Y C]);

        rss = sqrt(sum(abs(L).^2, 3));
        rss = single(rss) + eps('single');

        % Body mask from low-res RSS
        p95 = prctile(rss(:), 95);
        thr = max(eps('single'), maskFrac * p95);
        m = rss > thr;

        % Clean up mask a bit (fill holes + light dilation)
        m = imfill(m, 'holes');
        m = imdilate(m, strel('disk', 2));

        % Smooth RSS -> bias estimate
        rss_s = imgaussfilt(rss, smoothSigma, 'FilterSize', filtSz);

        % Normalise median in mask to 1
        medv = median(rss_s(m));
        if ~isfinite(medv) || medv <= 0, medv = 1; end
        b = rss_s / medv;

        bias(:,:,s)  = b;
        mask3(:,:,s) = m;
    end
end
