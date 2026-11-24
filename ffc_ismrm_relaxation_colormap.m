function [xDisplay, cmap] = ffc_ismrm_relaxation_colormap(mapType, x, loLev, upLev, useLogRemap)
%FFC_ISMRM_RELAXATION_COLORMAP ISMRM-recommended T1/R1 colour display.
%   [XDISPLAY, CMAP] = FFC_ISMRM_RELAXATION_COLORMAP(TYPE,X,LOW,HIGH,USELOG)
%   uses the Lipari scientific colour map recommended for T1/R1 relaxometry.
%   When USELOG is true (default), the colour table is remapped with the
%   logarithm-like mapping described by Fuderer et al., MRM 2025;93:490-506.
%   Image values remain in their native physical units; the logarithm is NOT
%   applied to the T1/R1 values themselves.
%
%   TYPE is 'T1' or 'R1'.  For R1 the Lipari table is reversed, matching the
%   MATLAB reference implementation distributed with the recommendation.
%   Zero/non-finite values are treated as invalid and map to black.

    if nargin < 5 || isempty(useLogRemap), useLogRemap = true; end
    mapType = upper(char(string(mapType)));
    if ~ismember(mapType, {'T1','R1'})
        error('ffc_ismrm_relaxation_colormap:BadMapType', ...
            'mapType must be ''T1'' or ''R1''.');
    end
    validateattributes(loLev, {'numeric'}, {'scalar','real','finite'});
    validateattributes(upLev, {'numeric'}, {'scalar','real','finite','>',loLev});
    if upLev <= 0
        error('ffc_ismrm_relaxation_colormap:BadRange', ...
            'Upper display level must be positive.');
    end

    persistent lipari
    if isempty(lipari)
        root = fileparts(mfilename('fullpath'));
        fn = fullfile(root,'resources','ismrm_colormaps','lipari.csv');
        lipari = dlmread(fn, ' ');
        if size(lipari,2) ~= 3 || size(lipari,1) < 2
            error('ffc_ismrm_relaxation_colormap:BadPalette', ...
                'Invalid Lipari colour-map resource: %s', fn);
        end
    end

    baseMap = lipari;
    if strcmp(mapType,'R1')
        baseMap = flipud(baseMap);
    end
    baseMap(1,:) = 0; % ISMRM invalid-value colour

    if useLogRemap
        cmap = localLogRemap(baseMap, loLev, upLev);
    else
        cmap = baseMap;
    end

    % Prepare display values so invalid zero/non-finite values occupy the
    % first (black) colour while valid sub-window values use the darkest
    % valid colour rather than being mistaken for invalid pixels.
    n = size(baseMap,1);
    delta = (upLev - loLev) / n;
    invalidValue = loLev - delta;
    minValidValue = loLev + 1.5*delta;

    xDisplay = double(x);
    invalid = ~isfinite(xDisplay) | xDisplay <= 0;
    lowValid = ~invalid & xDisplay < (loLev + delta);
    xDisplay(invalid) = invalidValue;
    xDisplay(lowValid) = minValidValue;
end

function out = localLogRemap(in, loLev, upLev)
% Log-like colour-table remapping used by the ISMRM reference implementation.
    n = size(in,1);
    a = exp(-1) * upLev;
    m = max(a, loLev);
    b = (1/n) + double(a >= loLev) * ((a-loLev)/(2*a-loLev));
    b = b + 1e-7;

    out = zeros(size(in));
    out(1,:) = in(1,:);
    logPortion = 1/(log(m)-log(upLev));

    for g = 2:n
        xv = g*(upLev-loLev)/n + loLev;
        if xv > m
            f = n*((log(m)-log(xv))*logPortion*(1-b)+b);
        elseif (loLev < a) && (xv > loLev)
            f = n*((xv-loLev)/(a-loLev)*(b-1/n)) + 1;
        else
            f = 1;
        end
        out(g,:) = in(min(n,1+floor(f)),:);
    end
end
