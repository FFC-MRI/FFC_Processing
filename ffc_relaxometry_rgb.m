function rgb = ffc_relaxometry_rgb(mapType, imageData, loLev, upLev, useLogRemap)
%FFC_RELAXOMETRY_RGB Render a T1/R1 map with the same ISMRM colours as GUI.
    if nargin < 5, useLogRemap = true; end
    [x, cmap] = ffc_ismrm_relaxation_colormap(mapType, imageData, ...
        loLev, upLev, useLogRemap);
    n = size(cmap,1);
    idx = 1 + floor((x-loLev) ./ (upLev-loLev) .* (n-1));
    idx(x < loLev | ~isfinite(x)) = 1;
    idx = max(1,min(n,idx));
    rgb = ind2rgb(double(idx),cmap);
end
