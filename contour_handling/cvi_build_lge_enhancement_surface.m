function [surface,details] = cvi_build_lge_enhancement_surface(overlays,seriesContours,dicomInfo)
%CVI_BUILD_LGE_ENHANCEMENT_SURFACE Build a 3-D mesh of the saved cvi42 LGE overlay.
%
% cvi42's Tissue Characterization module stores its displayed enhancement
% mask in WorkImageOverlays.  Fully classified pixels are kept in
% PixelCoordinates and the partial-pixel edge in SubpixelCoordinates.  This
% function reconstructs that sparse mask at the saved sub-pixel resolution,
% builds an isosurface, and converts the vertices into DICOM patient space.
%
% The WorkImageOverlay Slice index is tied to cvi42's image index.  For LGE
% workspaces it can be linked unambiguously to the same-frame
% saEnhancementReferenceMyoContour (the FWHM enhanced-area reference ROI).

    surface = struct('Type','','PhaseIndex',NaN,'Vertices',[],'Faces',[], ...
        'SliceCount',0,'SliceIndices',[]);
    details = struct('available',false,'overlayCount',0,'matchedOverlayCount',0, ...
        'unmatchedOverlayCount',0,'subpixelResolution',NaN);
    if isempty(overlays) || isempty(seriesContours), return; end

    anchors = seriesContours(strcmp({seriesContours.Type},'LGEEnhancementReference'));
    if isempty(anchors), return; end

    mapped = struct('overlay',{},'geometry',{},'sliceIndex',{},'frameIndex0',{});
    for k = 1:numel(overlays)
        ov = overlays(k);
        if ~isfinite(ov.SliceIndex), continue; end
        % In cvi42 Tissue Characterization workspaces examined so far the
        % WorkImageOverlay Slice value equals the frame suffix in the image
        % key (e.g. Slice 6 <-> ...multi6).  Anchor against an actual saved
        % LGE ROI rather than assuming a DICOM frame ordering globally.
        a = find([anchors.FrameIndex0] == ov.SliceIndex,1);
        if isempty(a), continue; end
        frameNumber = anchors(a).FrameIndex0 + 1;
        try
            g = cvi_dicom_frame_geometry(dicomInfo,frameNumber);
        catch
            continue;
        end
        item.overlay = ov;
        item.geometry = g;
        item.sliceIndex = anchors(a).SliceIndex;
        item.frameIndex0 = anchors(a).FrameIndex0;
        mapped(end+1,1) = item; %#ok<AGROW>
    end

    details.overlayCount = numel(overlays);
    details.matchedOverlayCount = numel(mapped);
    details.unmatchedOverlayCount = details.overlayCount - details.matchedOverlayCount;
    if numel(mapped) < 2, return; end

    % Sort by actual physical slice location, not frame number.  This makes
    % the result independent of whether the source stack is apex->base or
    % base->apex.
    normal = mapped(1).geometry.normalDirection(:);
    normal = normal ./ norm(normal);
    locations = zeros(numel(mapped),1);
    for k = 1:numel(mapped)
        locations(k) = dot(mapped(k).geometry.imagePositionPatient(:),normal);
    end
    [locations,order] = sort(locations);
    mapped = mapped(order);

    resolutions = arrayfun(@(x) double(x.overlay.SubpixelResolution),mapped);
    resolutions = resolutions(isfinite(resolutions) & resolutions >= 1);
    if isempty(resolutions), resolution = 1; else, resolution = round(max(resolutions)); end
    resolution = max(1,min(16,resolution));
    details.subpixelResolution = resolution;

    % Convert all positive overlay samples to one common integer sub-pixel
    % grid. PixelCoordinates are complete native pixels; expand each of them
    % to resolution x resolution cells, then merge the explicitly saved edge
    % subpixels. This reproduces the displayed cvi42 mask much more closely
    % than using the FWHM reference ROI itself.
    subPoints = cell(numel(mapped),1);
    allX = []; allY = [];
    for k = 1:numel(mapped)
        ov = mapped(k).overlay;
        r0 = max(1,round(double(ov.SubpixelResolution)));
        p = double(ov.PixelCoordinates(ov.PixelValues > 0,:));
        sp = double(ov.SubpixelCoordinates(ov.SubpixelValues > 0,:));

        expanded = zeros(size(p,1)*resolution*resolution,2);
        row = 0;
        if ~isempty(p)
            for q = 1:size(p,1)
                baseX = round(p(q,1)*resolution);
                baseY = round(p(q,2)*resolution);
                [dx,dy] = meshgrid(0:resolution-1,0:resolution-1);
                n = numel(dx);
                expanded(row+(1:n),:) = [baseX+dx(:),baseY+dy(:)];
                row = row+n;
            end
        end
        expanded = expanded(1:row,:);
        if ~isempty(sp)
            if r0 ~= resolution
                sp = round(sp .* (resolution/r0));
            else
                sp = round(sp);
            end
        end
        pts = unique([expanded;sp],'rows');
        subPoints{k} = pts;
        if ~isempty(pts)
            allX = [allX;pts(:,1)]; %#ok<AGROW>
            allY = [allY;pts(:,2)]; %#ok<AGROW>
        end
    end
    if isempty(allX), return; end

    padXY = 2;
    minX = floor(min(allX))-padXY; maxX = ceil(max(allX))+padXY;
    minY = floor(min(allY))-padXY; maxY = ceil(max(allY))+padXY;
    nx = maxX-minX+1; ny = maxY-minY+1;
    nz = numel(mapped)+2; % zero layers close the mesh at stack ends
    mask = false(ny,nx,nz);
    for k = 1:numel(mapped)
        pts = subPoints{k};
        if isempty(pts), continue; end
        x = round(pts(:,1)-minX+1);
        y = round(pts(:,2)-minY+1);
        keep = x>=1 & x<=nx & y>=1 & y<=ny;
        idx = sub2ind([ny nx nz],y(keep),x(keep),repmat(k+1,nnz(keep),1));
        mask(idx) = true;
    end

    if nnz(mask) < 8, return; end
    fv = isosurface(double(mask),0.5);
    if isempty(fv.vertices) || isempty(fv.faces), return; end
    % The saved sub-pixel mask can create a dense mesh. Reduce it before
    % repeated plane intersection so the live registration sliders remain
    % responsive on older MATLAB releases. Keep substantially more detail
    % than the 128-ish FCI display matrix can resolve.
    if size(fv.faces,1) > 30000
        try
            fv = reducepatch(fv,30000/size(fv.faces,1));
        catch
        end
    end

    % Extrapolate one zero layer at each end solely to close the binary
    % volume.  Patient position within the real stack is linearly
    % interpolated from the actual per-frame ImagePositionPatient values.
    if numel(locations) > 1
        d0 = median(diff(locations));
        if ~isfinite(d0) || d0 == 0, d0 = 1; end
    else
        d0 = 1;
    end
    zLocations = [locations(1)-d0; locations(:); locations(end)+d0];
    origins = zeros(nz,3);
    for k = 1:numel(mapped)
        origins(k+1,:) = mapped(k).geometry.imagePositionPatient(:).';
    end
    origins(1,:) = origins(2,:) - d0*normal.';
    origins(end,:) = origins(end-1,:) + d0*normal.';

    v = double(fv.vertices); % columns are x(column), y(row), z(slice)
    subX = minX + (v(:,1)-1);
    subY = minY + (v(:,2)-1);
    pixelX = subX ./ resolution;
    pixelY = subY ./ resolution;

    zq = v(:,3);
    originAtZ = zeros(size(v,1),3);
    gridZ = (1:nz).';
    for dim = 1:3
        originAtZ(:,dim) = interp1(gridZ,origins(:,dim),zq,'linear','extrap');
    end

    % The LGE stack is expected to be parallel. Average tiny per-frame
    % numerical differences in orientation and spacing.
    rowDirs = zeros(numel(mapped),3); colDirs = rowDirs;
    rowSpacing = zeros(numel(mapped),1); colSpacing = rowSpacing;
    for k = 1:numel(mapped)
        rowDirs(k,:) = mapped(k).geometry.rowDirection(:).';
        colDirs(k,:) = mapped(k).geometry.columnDirection(:).';
        rowSpacing(k) = mapped(k).geometry.rowSpacingMm;
        colSpacing(k) = mapped(k).geometry.columnSpacingMm;
    end
    rowDir = mean(rowDirs,1); rowDir = rowDir./norm(rowDir);
    colDir = mean(colDirs,1); colDir = colDir-rowDir*dot(rowDir,colDir); colDir = colDir./norm(colDir);
    dr = mean(rowSpacing,'omitnan'); dc = mean(colSpacing,'omitnan');

    patient = originAtZ + (pixelX*dc).*rowDir + (pixelY*dr).*colDir;
    surface.Type = 'LGEEnhancement';
    surface.PhaseIndex = 1;
    surface.Vertices = patient;
    surface.Faces = fv.faces;
    surface.SliceCount = numel(mapped);
    surface.SliceIndices = [mapped.sliceIndex];
    details.available = true;
end
