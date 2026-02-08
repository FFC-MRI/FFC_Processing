function img = ifftnc_spatial(k, nSpatial)
%IFFTNC_SPATIAL Centered ifftn over the first nSpatial dimensions, leaving
%other dimensions (slice/time/field/coil) intact.

    if nargin < 2, nSpatial = 3; end

    sz = size(k);
    sz(end+1:nSpatial) = 1;

    % Build fftshift/ifftshift vectors for the spatial dims only
    img = k;

    % ifftshift spatial dims
    for d = 1:nSpatial
        img = ifftshift(img, d);
    end

    % ifftn over spatial dims (MATLAB ifftn operates over all dims unless told;
    % so we reshape spatial dims into one block, do ifftn, then reshape back)
    % Easiest robust approach: apply ifft along each spatial dim:
    for d = 1:nSpatial
        img = ifft(img, [], d);
    end

    % fftshift spatial dims
    for d = 1:nSpatial
        img = fftshift(img, d);
    end
end
