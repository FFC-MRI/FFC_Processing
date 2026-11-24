function [mask, details] = cvi_series_plane_binary_mask(series, phaseIndex, ...
    sourceCenter, view, registration, imageSize, sourceName)
%CVI_SERIES_PLANE_BINARY_MASK Rasterise a registered CVI contour into an FCI plane.
%
% [MASK,DETAILS] = CVI_SERIES_PLANE_BINARY_MASK(SERIES,PHASEINDEX,
% SOURCECENTER,VIEW,REGISTRATION,IMAGESIZE,SOURCENAME) returns a logical
% mask in the currently DISPLAYED FCI plane.  The same registration,
% in-plane rotation and display-mirror corrections used for the contour
% overlay are applied before rasterisation.
%
% SOURCENAME may be:
%   'LGEEnhancement'           saved cvi42 LGE enhancement overlay
%   'LVMyocardium'             LV epicardium minus LV endocardium
%   'LVCavity'                 LV endocardial interior
%   'RVCavity'                 RV endocardial interior
%   'LGEEnhancementReference'  enhanced-area/FWHM reference ROI
%   'LGERemoteReference'       remote-myocardium reference ROI
%
% The returned mask is deliberately in display coordinates so callers can
% use ffc_restore_native_plane_mask to obtain the native reconstruction ROI
% used by ImageReconCore.T1dispersion.

    if nargin < 7 || isempty(sourceName)
        sourceName = 'LGEEnhancement';
    end
    sourceName = char(string(sourceName));
    if numel(imageSize) < 2
        error('cvi_series_plane_binary_mask:InvalidImageSize', ...
            'imageSize must contain [rows columns].');
    end
    rows = round(double(imageSize(1)));
    cols = round(double(imageSize(2)));
    if rows < 1 || cols < 1
        error('cvi_series_plane_binary_mask:InvalidImageSize', ...
            'The current FCI image has an invalid matrix size.');
    end

    details = struct('source',sourceName,'sourceLabel',maskSourceLabel(sourceName), ...
        'pixelCount',0,'areaMm2',0,'componentCount',0,'boundarySegmentCount',0);

    switch sourceName
        case 'LVMyocardium'
            epi = surfaceMask('LVEpi',true);
            endo = surfaceMask('LVEndo',true);
            if ~any(epi(:))
                error('cvi_series_plane_binary_mask:NoLVEpiIntersection', ...
                    'The LV epicardial surface does not intersect the current FCI slice.');
            end
            if ~any(endo(:))
                error('cvi_series_plane_binary_mask:NoLVEndoIntersection', ...
                    'The LV endocardial surface does not intersect the current FCI slice.');
            end
            mask = epi & ~endo;
        case 'LVCavity'
            mask = surfaceMask('LVEndo',true);
        case 'RVCavity'
            mask = surfaceMask('RVEndo',true);
        case 'LGEEnhancement'
            mask = surfaceMask('LGEEnhancement',false);
        case 'LGEEnhancementReference'
            mask = surfaceMask('LGEEnhancementReference',true);
        case 'LGERemoteReference'
            mask = surfaceMask('LGERemoteReference',true);
        otherwise
            error('cvi_series_plane_binary_mask:UnknownSource', ...
                'Unknown CVI binary-mask source "%s".',sourceName);
    end

    mask = logical(mask);
    if ~any(mask(:))
        error('cvi_series_plane_binary_mask:EmptyMask', ...
            '%s does not intersect the current FCI slice.',details.sourceLabel);
    end

    details.pixelCount = nnz(mask);
    pixelX = double(view.fovX) ./ cols;
    pixelY = double(view.fovY) ./ rows;
    details.areaMm2 = details.pixelCount .* pixelX .* pixelY;
    try
        cc = bwconncomp(mask);
        details.componentCount = cc.NumObjects;
    catch
        details.componentCount = NaN;
    end

    function out = surfaceMask(typeName,phaseSpecific)
        if ~isfield(series,'Surfaces') || isempty(series.Surfaces)
            error('cvi_series_plane_binary_mask:NoSurfaces', ...
                'The selected CVI series contains no reconstructed contour surfaces.');
        end
        match = strcmp({series.Surfaces.Type},typeName);
        if phaseSpecific
            match = match & [series.Surfaces.PhaseIndex] == phaseIndex;
        end
        indices = find(match);
        if isempty(indices)
            error('cvi_series_plane_binary_mask:SurfaceUnavailable', ...
                '%s is not available for the selected clinical series/phase.', ...
                maskSourceLabel(typeName));
        end
        out = false(rows,cols);
        for si = indices(:).'
            [x,y] = cvi_surface_plane_intersection(series.Surfaces(si), ...
                view,registration,sourceCenter);
            if isempty(x), continue; end
            details.boundarySegmentCount = details.boundarySegmentCount + ...
                countSegments(x,y);
            out = out | segmentsToFilledMask(x,y,rows,cols,view.fovX,view.fovY);
        end
    end
end


function mask = segmentsToFilledMask(x,y,rows,cols,fovX,fovY)
% Rasterise unordered triangle/plane intersection segments.  The segments
% from a closed contour surface meet at their endpoints; over-sampling each
% segment onto the image grid gives a closed digital boundary which imfill
% can convert into the binary interior.  A small morphological close is used
% only as a fallback for sub-pixel endpoint rounding gaps.

    x = double(x(:)); y = double(y(:));
    edge = false(rows,cols);
    dx = double(fovX) ./ cols;
    dy = double(fovY) ./ rows;

    finite = isfinite(x) & isfinite(y);
    starts = find(finite & [true; ~finite(1:end-1)]);
    stops  = find(finite & [~finite(2:end); true]);
    for k = 1:min(numel(starts),numel(stops))
        idx = starts(k):stops(k);
        if numel(idx) < 2, continue; end
        % cvi_surface_plane_intersection currently emits two points per
        % triangle segment.  Supporting longer runs costs nothing and keeps
        % this helper tolerant of future ordered-polyline output.
        for q = 1:numel(idx)-1
            x1=x(idx(q)); x2=x(idx(q+1));
            y1=y(idx(q)); y2=y(idx(q+1));
            c1=(x1 + double(fovX)/2)./dx + 0.5;
            c2=(x2 + double(fovX)/2)./dx + 0.5;
            r1=(y1 + double(fovY)/2)./dy + 0.5;
            r2=(y2 + double(fovY)/2)./dy + 0.5;
            n = max(2,ceil(max(abs(c2-c1),abs(r2-r1))*4)+1);
            cc = round(linspace(c1,c2,n));
            rr = round(linspace(r1,r2,n));
            keep = cc>=1 & cc<=cols & rr>=1 & rr<=rows;
            if any(keep)
                edge(sub2ind([rows cols],rr(keep),cc(keep))) = true;
            end
        end
    end

    if ~any(edge(:))
        mask = edge;
        return;
    end

    mask = imfill(edge,'holes');
    % If nothing appreciably larger than the one-pixel boundary was filled,
    % endpoint quantisation probably left a one-pixel gap.  Close only in
    % that case so normal masks retain the original CVI boundary as closely
    % as possible.
    if nnz(mask) <= nnz(edge) + max(4,round(0.02*nnz(edge)))
        try
            closed = imclose(edge,strel('square',3));
            candidate = imfill(closed,'holes');
            if nnz(candidate) > nnz(mask)
                mask = candidate;
            end
        catch
            % Image Processing Toolbox is already required elsewhere in the
            % FFC suite, but retain a harmless fallback for older releases.
        end
    end
end


function n = countSegments(x,y)
    finite = isfinite(x(:)) & isfinite(y(:));
    n = nnz(finite & [true; ~finite(1:end-1)]);
end


function label = maskSourceLabel(sourceName)
    switch char(string(sourceName))
        case 'LGEEnhancement',          label = 'LGE enhancement';
        case 'LVMyocardium',            label = 'LV myocardium';
        case 'LVCavity',                label = 'LV blood pool';
        case 'RVCavity',                label = 'RV blood pool';
        case 'LGEEnhancementReference', label = 'Enhanced-area reference ROI';
        case 'LGERemoteReference',      label = 'Remote myocardium ROI';
        case 'LVEpi',                   label = 'LV epicardium';
        case 'LVEndo',                  label = 'LV endocardium';
        case 'RVEndo',                  label = 'RV endocardium';
        otherwise,                      label = char(string(sourceName));
    end
end
