function geometry = ffc_scan_spatial_geometry(scan)
%FFC_SCAN_SPATIAL_GEOMETRY Patient-space geometry for imaging or SVS data.
%
% Array dimensions are read, phase and slice/partition.  IMAGE_ORIENTATION_
% SUBJECT supplies the corresponding read and phase unit vectors; their
% cross-product supplies the prescription normal. Returned patient coordinates
% are DICOM LPS in mm. Cameleon scanner coordinates are converted once here.
% TEST25 Y/Z conversion is based on the supplied HFS axial/SAX datasets.

    orientation = memberValue(scan, 'orientation', []);
    if isempty(orientation)
        orientation = parameterValue(scan, 'IMAGE_ORIENTATION_SUBJECT', []);
    end

    if isnumeric(orientation) && numel(orientation) == 6 && ...
            all(isfinite(double(orientation(:))))
        readDirection = double(orientation(1:3));
        prescriptionPhaseDirection = double(orientation(4:6));
        readDirection = readDirection(:) ./ norm(readDirection);
        prescriptionPhaseDirection = prescriptionPhaseDirection(:) - ...
            readDirection .* dot(readDirection, prescriptionPhaseDirection(:));
        if norm(prescriptionPhaseDirection) < eps
            prescriptionPhaseDirection = [0; -1; 0];
        end
        prescriptionPhaseDirection = prescriptionPhaseDirection ./ ...
            norm(prescriptionPhaseDirection);
    else
        readDirection = [1; 0; 0];
        prescriptionPhaseDirection = [0; -1; 0];
    end
    prescriptionSliceDirection = cross(readDirection, ...
        prescriptionPhaseDirection);
    prescriptionSliceDirection = prescriptionSliceDirection ./ ...
        norm(prescriptionSliceDirection);

    fovRead = positiveScalar(memberValue(scan, 'fov', []), NaN);
    fovPhase = positiveScalar(memberValue(scan, 'fov_phase', []), fovRead);
    fovSlice = positiveScalar(memberValue(scan, 'fov_3d', []), NaN);
    if ~isfinite(fovRead)
        fovRead = 1000 .* positiveScalar( ...
            parameterValue(scan, 'FIELD_OF_VIEW', []), 1);
    end
    if ~isfinite(fovPhase)
        fovPhase = 1000 .* positiveScalar( ...
            parameterValue(scan, 'FIELD_OF_VIEW_PHASE', []), fovRead ./ 1000);
    end
    if ~isfinite(fovSlice)
        fovSlice = 1000 .* positiveScalar( ...
            parameterValue(scan, 'FIELD_OF_VIEW_3D', []), NaN);
    end
    if ~isfinite(fovSlice)
        thickness = positiveScalar(memberValue(scan, 'thk', []), 1);
        fovSlice = thickness .* max(positiveScalar( ...
            memberValue(scan, 'slices', []), 1), 1);
    end

    offsetRead = 1000 .* finiteScalar( ...
        parameterValue(scan, 'OFF_CENTER_FIELD_OF_VIEW_1D', 0), 0);
    offsetPhase = 1000 .* finiteScalar( ...
        parameterValue(scan, 'OFF_CENTER_FIELD_OF_VIEW_2D', 0), 0);
    offsetSlice = 1000 .* finiteScalar( ...
        parameterValue(scan, 'OFF_CENTER_FIELD_OF_VIEW_3D', 0), 0);

    % Spinlab stores the prescription centre twice: once in the local
    % read/phase/slice frame and once directly in scanner X/Y/Z.  The latter
    % is the position shown by the planning UI and avoids compounding any
    % row/column sign convention ambiguity in IMAGE_ORIENTATION_SUBJECT.
    % Prefer it when the complete triplet is present; retain the local-frame
    % calculation as a fallback for older datasets.
    prescriptionBasis = [readDirection, prescriptionPhaseDirection, ...
        prescriptionSliceDirection];
    localFrameCenterMm = prescriptionBasis * ...
        [offsetRead; offsetPhase; offsetSlice];
    scannerCenterM = [ ...
        finiteParameterOrNaN(scan, 'OFF_CENTER_FIELD_OF_VIEW_X'); ...
        finiteParameterOrNaN(scan, 'OFF_CENTER_FIELD_OF_VIEW_Y'); ...
        finiteParameterOrNaN(scan, 'OFF_CENTER_FIELD_OF_VIEW_Z')];
    if all(isfinite(scannerCenterM))
        centerPatientMm = 1000 .* scannerCenterM;
        centerSource = 'OFF_CENTER_FIELD_OF_VIEW_X/Y/Z';
    else
        centerPatientMm = localFrameCenterMm;
        centerSource = 'oriented read/phase/slice offsets';
    end

    imageData = memberValue(scan, 'magimage', []);
    if isempty(imageData)
        matrixSize = [positiveScalar(memberValue(scan, 'samples', []), 1), ...
            positiveScalar(memberValue(scan, 'views', []), 1), ...
            positiveScalar(memberValue(scan, 'slices', []), 1)];
    else
        matrixSize = [size(imageData, 1), size(imageData, 2), ...
            size(imageData, 3)];
    end
    matrixSize = max(1, round(double(matrixSize)));

    twoDimensional = memberValue(scan, 'TwoDimensional', 1);
    isSVS = logicalScalar(memberValue(scan, 'is_svs', false), false);

    % Only Cameleon metadata identifies the scanner convention reviewed in
    % TEST25. Do not silently apply this conversion to unrelated file formats.
    parameters = memberValue(scan, 'param', struct());
    isCameleon = isstruct(parameters) && ...
        isfield(parameters, 'IMAGE_ORIENTATION_SUBJECT');
    scannerToPatient = eye(3);
    if isCameleon
        scannerToPatient = diag([1 -1 -1]);
        % The empirical TEST24 export maps native [read phase slice] to
        % DICOM by reversing both in-plane indices and retaining slice order.
        % Use the SAME mapping for display/contours, not the previous phase-
        % and-slice reversal. 3-D ordering still needs scanner validation.
        arrayBasis = prescriptionBasis * diag([-1 -1 1]);
    else
        % Preserve the existing array mapping for legacy input formats.
        arrayBasis = prescriptionBasis * diag([1 -1 -1]);
    end
    geometry = struct();
    geometry.dimensionBasisPatient = scannerToPatient * arrayBasis;
    geometry.prescriptionBasisPatient = scannerToPatient * prescriptionBasis;
    geometry.scannerToPatient = scannerToPatient;
    geometry.coordinateConvention = 'DICOM LPS';
    geometry.geometryVersion = 25;
    geometry.centerPatientMm = scannerToPatient * centerPatientMm;
    geometry.centerSource = centerSource;
    geometry.localFrameCenterPatientMm = scannerToPatient * localFrameCenterMm;
    geometry.centerConsistencyErrorMm = norm(centerPatientMm - localFrameCenterMm);
    geometry.fovMm = [fovRead, fovPhase, fovSlice];
    geometry.matrixSize = matrixSize;
    geometry.spacingMm = geometry.fovMm ./ matrixSize;
    geometry.is3D = ~isSVS && finiteScalar(twoDimensional, 1) ~= 1 && ...
        matrixSize(3) > 1;
    if ~isSVS && finiteScalar(twoDimensional, 1) == 1
        thicknessMm = 1000 * positiveScalar( ...
            parameterValue(scan, 'SLICE_THICKNESS', []), NaN);
        if ~isfinite(thicknessMm)
            thicknessMm = positiveScalar(memberValue(scan, 'thk', []), NaN);
        end
        gapMm = 1000 * finiteScalar( ...
            parameterValue(scan, 'SPACING_BETWEEN_SLICE', 0), 0);
        if gapMm < 0
            error('ffc_scan_spatial_geometry:NegativeGap', ...
                'SPACING_BETWEEN_SLICE is a physical gap and must be nonnegative.');
        end
        if isfinite(thicknessMm)
            geometry.spacingMm(3) = thicknessMm + gapMm;
        end
    end
end


function value = finiteParameterOrNaN(scan, name)
    value = finiteScalar(parameterValue(scan, name, NaN), NaN);
end


function value = memberValue(scan, name, defaultValue)
    value = defaultValue;
    try
        if (isobject(scan) && isprop(scan, name)) || ...
                (isstruct(scan) && isfield(scan, name))
            candidate = scan.(name);
            if ~isempty(candidate)
                value = candidate;
            end
        end
    catch
    end
end


function value = parameterValue(scan, name, defaultValue)
    parameters = memberValue(scan, 'param', struct());
    value = defaultValue;
    if isstruct(parameters) && isfield(parameters, name) && ...
            ~isempty(parameters.(name))
        value = parameters.(name);
    end
end


function value = positiveScalar(candidate, defaultValue)
    value = finiteScalar(candidate, defaultValue);
    if ~isfinite(value) || value <= 0
        value = defaultValue;
    end
end


function value = finiteScalar(candidate, defaultValue)
    value = defaultValue;
    if isnumeric(candidate) || islogical(candidate)
        candidate = double(candidate(:));
        candidate = candidate(isfinite(candidate));
        if ~isempty(candidate)
            value = candidate(1);
        end
    end
end


function value = logicalScalar(candidate, defaultValue)
    value = defaultValue;
    if (isnumeric(candidate) || islogical(candidate)) && ...
            isscalar(candidate) && isfinite(double(candidate))
        value = logical(candidate);
    end
end
