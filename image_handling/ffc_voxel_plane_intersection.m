function [x, y, details] = ffc_voxel_plane_intersection(svsScan, view)
%FFC_VOXEL_PLANE_INTERSECTION Intersect an SVS cuboid with a displayed plane.

    svsGeometry = ffc_scan_spatial_geometry(svsScan);
    voxelSize = [];
    try
        if isobject(svsScan) && isprop(svsScan, 'voxel_size_mm')
            voxelSize = svsScan.voxel_size_mm;
        elseif isstruct(svsScan) && isfield(svsScan, 'voxel_size_mm')
            voxelSize = svsScan.voxel_size_mm;
        end
    catch
    end
    if isempty(voxelSize) || numel(voxelSize) ~= 3 || ...
            any(~isfinite(double(voxelSize))) || any(double(voxelSize) <= 0)
        voxelSize = svsGeometry.fovMm;
    end
    voxelSize = double(voxelSize(:));

    signs = [-1 -1 -1; 1 -1 -1; -1 1 -1; 1 1 -1; ...
             -1 -1 1; 1 -1 1; -1 1 1; 1 1 1].';
    voxelBasis = svsGeometry.dimensionBasisPatient;
    if isfield(svsGeometry, 'prescriptionBasisPatient')
        voxelBasis = svsGeometry.prescriptionBasisPatient;
    end
    corners = svsGeometry.centerPatientMm + ...
        voxelBasis * ...
        (0.5 .* voxelSize .* signs);
    edges = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; 3 7; ...
        4 8; 5 6; 5 7; 6 8; 7 8];

    planeNormal = double(view.basisNormalPatient(:));
    planePoint = double(view.slicePointPatientMm(:));
    signedDistance = planeNormal.' * (corners - planePoint);
    tolerance = max(voxelSize) .* 1e-9 + eps;
    points = corners(:, abs(signedDistance) <= tolerance);

    for edgeIndex = 1:size(edges, 1)
        firstIndex = edges(edgeIndex, 1);
        secondIndex = edges(edgeIndex, 2);
        firstDistance = signedDistance(firstIndex);
        secondDistance = signedDistance(secondIndex);
        if firstDistance .* secondDistance < 0
            fraction = firstDistance ./ (firstDistance - secondDistance);
            points(:, end + 1) = corners(:, firstIndex) + fraction .* ...
                (corners(:, secondIndex) - corners(:, firstIndex)); %#ok<AGROW>
        end
    end

    points = uniqueColumns(points, tolerance);
    x = [];
    y = [];
    details = struct('intersects', false, 'patientPointsMm', points, ...
        'voxelCenterPatientMm', svsGeometry.centerPatientMm, ...
        'voxelSizeMm', voxelSize(:).');
    if size(points, 2) < 3
        return;
    end

    relative = points - view.geometry.centerPatientMm;
    xUnordered = view.basisXPatient.' * relative;
    yUnordered = view.basisYPatient.' * relative;
    hull = convhull(xUnordered(:), yUnordered(:));
    x = xUnordered(hull);
    y = yUnordered(hull);
    details.intersects = true;
end


function uniquePoints = uniqueColumns(points, tolerance)
    uniquePoints = zeros(3, 0);
    for pointIndex = 1:size(points, 2)
        candidate = points(:, pointIndex);
        if isempty(uniquePoints) || all(vecnorm(uniquePoints - candidate, 2, 1) > tolerance)
            uniquePoints(:, end + 1) = candidate; %#ok<AGROW>
        end
    end
end
