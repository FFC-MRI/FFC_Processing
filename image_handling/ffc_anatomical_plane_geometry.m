function view = ffc_anatomical_plane_geometry(scan, planeName, sliceIndex, geomT)
%FFC_ANATOMICAL_PLANE_GEOMETRY Map a volume onto a radiological screen.
%
% The returned basis vectors describe increasing CData columns (basisX) and
% increasing CData rows (basisY).  MATLAB imagesc displays rows from top to
% bottom, so the compass uses -basisY for its screen-up direction.  The
% acquisition dimensions are selected from IMAGE_ORIENTATION_SUBJECT and are
% then signed to the same patient convention for every 2D/3D acquisition.
% RS2D subject coordinates are +X=R, +Y=A and +Z=H, so the canonical
% radiological display has axial/coronal left=R, right=L and axial
% top=A, bottom=P.

    if nargin < 2 || isempty(planeName)
        planeName = 'Axial';
    end
    if nargin < 3 || isempty(sliceIndex)
        sliceIndex = 1;
    end
    if nargin < 4 || isempty(geomT)
        geomT = eye(3);
    end

    geometry = ffc_scan_spatial_geometry(scan);
    basis = geometry.dimensionBasisPatient;
    planeName = validatestring(char(string(planeName)), ...
        {'Axial', 'Coronal', 'Sagittal'});

    switch planeName
        case 'Axial'
            targetNormal = [0; 0; 1];
            targetHorizontal = [-1; 0; 0]; % screen right = patient L
            targetRowDown = [0; -1; 0];    % screen down = patient P
        case 'Coronal'
            targetNormal = [0; 1; 0];
            targetHorizontal = [-1; 0; 0]; % screen right = patient L
            targetRowDown = [0; 0; -1];    % screen down = patient F
        otherwise % Sagittal
            targetNormal = [1; 0; 0];
            targetHorizontal = [0; -1; 0]; % screen right = patient P
            targetRowDown = [0; 0; -1];    % screen down = patient F
    end

    if geometry.is3D
        [~, normalDimension] = max(abs(basis.' * targetNormal));
        remaining = setdiff(1:3, normalDimension, 'stable');
        directScore = abs(dot(basis(:, remaining(1)), targetHorizontal)) + ...
            abs(dot(basis(:, remaining(2)), targetRowDown));
        swappedScore = abs(dot(basis(:, remaining(2)), targetHorizontal)) + ...
            abs(dot(basis(:, remaining(1)), targetRowDown));
        if directScore >= swappedScore
            columnDimension = remaining(1);
            rowDimension = remaining(2);
        else
            columnDimension = remaining(2);
            rowDimension = remaining(1);
        end
    else
        % Multislice 2D data have an acquired plane and a slice label.  The
        % third array dimension is therefore the normal; only the two in-plane
        % dimensions may be exchanged to obtain the canonical screen.
        normalDimension = 3;
        remaining = [1 2];
        directScore = abs(dot(basis(:, remaining(2)), targetHorizontal)) + ...
            abs(dot(basis(:, remaining(1)), targetRowDown));
        swappedScore = abs(dot(basis(:, remaining(1)), targetHorizontal)) + ...
            abs(dot(basis(:, remaining(2)), targetRowDown));
        if directScore >= swappedScore
            columnDimension = remaining(2);
            rowDimension = remaining(1);
        else
            columnDimension = remaining(1);
            rowDimension = remaining(2);
        end
    end

    nSlices = geometry.matrixSize(normalDimension);
    sliceIndex = min(max(round(double(sliceIndex)), 1), nSlices);
    normalOffset = (sliceIndex - (nSlices + 1) ./ 2) .* ...
        geometry.spacingMm(normalDimension);

    nativeBasisX = basis(:, columnDimension);
    nativeBasisY = basis(:, rowDimension);
    basisNormal = basis(:, normalDimension);

    % Canonicalise the acquisition signs before applying user display
    % rotation/reflection.  The signed diagonal is also the exact lossless
    % image operation needed by ffc_extract_anatomical_plane.
    baseTransform = diag([signedAlignment(nativeBasisX, targetHorizontal), ...
        signedAlignment(nativeBasisY, targetRowDown)]);
    userTransform = validateInplaneTransform(geomT);
    effectiveTransform = userTransform * baseTransform;

    transformedBasis = [nativeBasisX, nativeBasisY] * effectiveTransform;
    fovXY = [geometry.fovMm(columnDimension); geometry.fovMm(rowDimension)];
    transformedFov = abs(effectiveTransform).' * fovXY;

    view = struct();
    view.plane = planeName;
    view.rowDimension = rowDimension;
    view.columnDimension = columnDimension;
    view.normalDimension = normalDimension;
    view.nSlices = nSlices;
    view.sliceIndex = sliceIndex;
    view.basisXPatient = transformedBasis(:, 1);
    view.basisYPatient = transformedBasis(:, 2);
    view.basisNormalPatient = basisNormal;
    view.fovX = transformedFov(1);
    view.fovY = transformedFov(2);
    view.slicePointPatientMm = geometry.centerPatientMm + ...
        normalOffset .* basisNormal;
    view.spacingNormalMm = geometry.spacingMm(normalDimension);
    view.geometry = geometry;
    view.geomT = geomT;
    view.baseTransform = baseTransform;
    view.effectiveTransform = effectiveTransform;
end


function signValue = signedAlignment(candidate, target)
    if dot(candidate, target) < 0
        signValue = -1;
    else
        signValue = 1;
    end
end


function transform = validateInplaneTransform(geomT)
    if ~isnumeric(geomT) || size(geomT, 1) < 2 || size(geomT, 2) < 2
        error('ffc_anatomical_plane_geometry:InvalidTransform', ...
            'geomT must contain a signed in-plane permutation.');
    end
    raw = double(geomT(1:2, 1:2));
    transform = round(raw);
    if any(~isfinite(raw(:))) || max(abs(raw(:) - transform(:))) > 1e-8 || ...
            any(sum(abs(transform), 1) ~= 1) || ...
            any(sum(abs(transform), 2) ~= 1)
        error('ffc_anatomical_plane_geometry:InvalidTransform', ...
            'geomT must contain a signed in-plane permutation.');
    end
end
