function [plane, view] = ffc_extract_anatomical_plane(volume, scan, planeName, sliceIndex, geomT)
%FFC_EXTRACT_ANATOMICAL_PLANE Extract and losslessly orient one volume plane.

    if nargin < 5 || isempty(geomT)
        geomT = eye(3);
    end
    view = ffc_anatomical_plane_geometry(scan, planeName, sliceIndex, geomT);

    volumeSize = size(volume);
    volumeSize(end + 1:3) = 1;
    volume = reshape(volume, volumeSize);
    permutation = [view.rowDimension, view.columnDimension, ...
        view.normalDimension, 4:ndims(volume)];
    % Select the thin native slab before permuting, avoiding a volume copy.
    indices = repmat({':'}, 1, max(3, ndims(volume)));
    indices{view.normalDimension} = min(view.sliceIndex, size(volume,view.normalDimension));
    permuted = permute(volume(indices{:}), permutation);
    plane = permuted(:, :, 1, :);
    plane = reshape(plane, size(permuted, 1), size(permuted, 2), []);
    if size(plane, 3) == 1
        plane = plane(:, :, 1);
    end
    % The geometry helper includes both the acquisition-to-radiological
    % canonicalisation and the user display transform in effectiveTransform.
    plane = apply_inplane_geometry(plane, view.effectiveTransform);
end
