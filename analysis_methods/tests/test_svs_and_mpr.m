function tests = test_svs_and_mpr
% Regression tests for SVS reduction, MPR geometry and voxel overlays.

    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.originalPath = path;
    addpath(projectRoot);
    addpath(fullfile(projectRoot, 'image_handling'));
    addpath(fullfile(projectRoot, 'kspace_handling'));
end


function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end


function testSVSSignalUsesMagnitudeBeforeSampleMean(testCase)
    data = complex(zeros(4, 1, 1, 2, 3, 1, 'single'));
    for timeIndex = 1:2
        for fieldIndex = 1:3
            scale = single(timeIndex + 10 .* fieldIndex);
            values = scale .* single([1 -1 1 -1]);
            data(:, 1, 1, timeIndex, fieldIndex, 1) = complex(values, 0);
        end
    end

    [complexSignal, perReceiver, magnitudeSignal, ...
        perReceiverMagnitude] = svs_signal_from_kspace(data);
    expectedComplex = mean(data, 1);
    expectedMagnitude = mean(abs(data), 1);

    verifyEqual(testCase, size(complexSignal), [1 1 1 2 3]);
    verifyEqual(testCase, perReceiver, expectedComplex);
    verifyEqual(testCase, complexSignal, ...
        expectedComplex(:, :, :, :, :, 1));
    verifyEqual(testCase, perReceiverMagnitude, expectedMagnitude);
    verifyEqual(testCase, magnitudeSignal, ...
        expectedMagnitude(:, :, :, :, :, 1));
    verifyGreaterThan(testCase, min(magnitudeSignal(:)), 0);
    verifyEqual(testCase, complexSignal, zeros(size(complexSignal), 'single'));
end


function testSVSConstructorToleratesMissingImagingFields(testCase)
    raw = complex(zeros(4, 1, 1, 6, 1, 'single'));
    for acquisition = 1:6
        raw(:, 1, 1, acquisition, 1) = complex(single(1:4), ...
            -single(1:4));
    end
    par = struct();
    par.cameleon = struct( ...
        'FIELD_OF_VIEW', 0.025, ...
        'FIELD_OF_VIEW_PHASE', 0.025, ...
        'FIELD_OF_VIEW_3D', 0.025, ...
        'SLICE_THICKNESS', 0.025, ...
        'IMAGE_ORIENTATION_SUBJECT', [1 0 0 0 -1 0], ...
        'FCI_LIST_OF_BEVO', [200 100], ...
        'FCI_NUMBER_OF_TEVO', 3, ...
        'FCI_LIST_OF_T1_ESTIMATES', [0.25 0.25], ...
        'FCI_T1_FACTOR_LOWER', 0.05, ...
        'FCI_T1_FACTOR_UPPER', 8, ...
        'FCI_T1_LOG_DIST', true);
    par.md1d = 4;
    par.md2d = 1;
    par.md3d = 1;
    par.md4d = 6;
    par.ncoils = 1;
    par.navg = 1;
    par.te = 0;
    par.fa = 90;
    par.serieName = 'Field-Cycling SVS';
    par.seqName = 'SVS';
    input = struct('data', raw, 'par', par);

    scan = ImageReconCore(input, '.dat', 1);
    scan = preprocessing(scan);
    scan = buildimages(scan);

    verifyTrue(testCase, scan.is_svs);
    verifyEqual(testCase, size(scan.compleximage), [1 1 1 3 2]);
    verifyEqual(testCase, scan.compleximage, ...
        reshape(mean(raw, 1), [1 1 1 3 2]), 'AbsTol', 1e-6);
    verifyEqual(testCase, scan.magimage, ...
        reshape(mean(abs(raw), 1), [1 1 1 3 2]), 'AbsTol', 1e-6);
    verifyEqual(testCase, scan.param.SVS_PROCESSED_MAGNITUDE_SIGNAL, ...
        reshape(mean(abs(raw), 1), [3 2]), 'AbsTol', 1e-6);
end


function testSimpleIrFitRecoversOffsetExponential(testCase)
    timesMs = [10 25 50 100 200 400 800].';
    expectedT1ms = 175;
    signedSignal = -240 .* exp(-timesMs ./ expectedT1ms) + 125;

    result = ffc_simple_ir_fit(timesMs, signedSignal);

    verifyEqual(testCase, result.T1ms, expectedT1ms, 'RelTol', 1e-5);
    verifyEqual(testCase, result.C, 125, 'RelTol', 1e-5);
    verifyEqual(testCase, result.A, -240, 'RelTol', 1e-5);
    verifyLessThan(testCase, result.sse, 1e-8);
end


function testIrPrefixIsSelectedFromAllFields(testCase)
    normalizedTimes = [0.05 0.15 0.35 0.60 0.90 1.5 3.0].';
    expectedT1ms = [120 170 230];
    timesMs = normalizedTimes .* expectedT1ms;
    signedSignals = -220 .* exp(-timesMs ./ expectedT1ms) + 120;
    magnitudeSignals = abs(signedSignals);

    [negativePrefixCount, scores, fits] = ...
        ffc_select_ir_prefix(timesMs, magnitudeSignals);

    verifyEqual(testCase, negativePrefixCount, 4);
    verifyEqual(testCase, numel(scores), size(timesMs, 1) - 1);
    verifyEqual(testCase, numel(fits), numel(expectedT1ms));
    recovered = cellfun(@(item) item.T1ms, fits);
    verifyEqual(testCase, recovered, expectedT1ms, 'RelTol', 1e-5);
end


function testScannerXYZPrescriptionCentreTakesPriority(testCase)
    scan = syntheticSVSScan();
    scan.param.OFF_CENTER_FIELD_OF_VIEW_1D = 0.101;
    scan.param.OFF_CENTER_FIELD_OF_VIEW_2D = 0.102;
    scan.param.OFF_CENTER_FIELD_OF_VIEW_3D = 0.103;
    scan.param.OFF_CENTER_FIELD_OF_VIEW_X = -0.02149;
    scan.param.OFF_CENTER_FIELD_OF_VIEW_Y = 0.02344;
    scan.param.OFF_CENTER_FIELD_OF_VIEW_Z = -0.02734;

    geometry = ffc_scan_spatial_geometry(scan);

    verifyEqual(testCase, geometry.centerPatientMm, ...
        [-21.49; 23.44; -27.34], 'AbsTol', 1e-10);
    verifyEqual(testCase, geometry.centerSource, ...
        'OFF_CENTER_FIELD_OF_VIEW_X/Y/Z');
end


function testOrthogonalPlanesUseExpectedDimensions(testCase)
    scan = syntheticImagingScan();
    volume = reshape(single(1:(4*6*8)), [4 6 8]);

    [axial, axialView] = ffc_extract_anatomical_plane( ...
        volume, scan, 'Axial', 4, eye(3));
    [coronal, coronalView] = ffc_extract_anatomical_plane( ...
        volume, scan, 'Coronal', 3, eye(3));
    [sagittal, sagittalView] = ffc_extract_anatomical_plane( ...
        volume, scan, 'Sagittal', 2, eye(3));

    verifyEqual(testCase, size(axial), [6 4]);
    verifyEqual(testCase, size(coronal), [8 4]);
    verifyEqual(testCase, size(sagittal), [8 6]);
    verifyEqual(testCase, axialView.normalDimension, 3);
    verifyEqual(testCase, coronalView.normalDimension, 2);
    verifyEqual(testCase, sagittalView.normalDimension, 1);
    verifyEqual(testCase, axialView.fovX, 40);
    verifyEqual(testCase, axialView.fovY, 60);
    verifyEqual(testCase, axialView.basisXPatient, [-1; 0; 0], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, axialView.basisYPatient, [0; -1; 0], ...
        'AbsTol', 1e-12);
end


function testDisplayRotationIsLosslessAndUpdatesCompassBasis(testCase)
    scan = syntheticImagingScan();
    volume = reshape(single(1:(4*6*8)), [4 6 8]);
    clockwise = [0 1 0; -1 0 0; 0 0 1];

    [native, nativeView] = ffc_extract_anatomical_plane( ...
        volume, scan, 'Axial', 4, eye(3));
    [rotated, rotatedView] = ffc_extract_anatomical_plane( ...
        volume, scan, 'Axial', 4, clockwise);

    verifyEqual(testCase, rotated, flip(native.', 2));
    verifyEqual(testCase, size(rotated), fliplr(size(native)));
    verifyEqual(testCase, rotatedView.basisXPatient, ...
        -nativeView.basisYPatient, 'AbsTol', 1e-12);
    verifyEqual(testCase, rotatedView.basisYPatient, ...
        nativeView.basisXPatient, 'AbsTol', 1e-12);
    verifyEqual(testCase, [rotatedView.fovX rotatedView.fovY], ...
        [nativeView.fovY nativeView.fovX], 'AbsTol', 1e-12);
end


function testRadiologicalSignsAreCanonicalForTwoDimensionalData(testCase)
    scan = syntheticImagingScan();
    scan.TwoDimensional = 1;
    scan.slices = 1;
    scan.magimage = zeros(4, 6, 1, 'single');

    view = ffc_anatomical_plane_geometry(scan, 'Axial', 1, eye(3));
    verifyEqual(testCase, view.basisXPatient, [-1; 0; 0], ...
        'AbsTol', 1e-12);
    verifyEqual(testCase, view.basisYPatient, [0; -1; 0], ...
        'AbsTol', 1e-12);
end


function testProvidedSVSIsPosteriorToSpinEchoCentre(testCase)
    imageScan = syntheticImagingScan();
    imageScan.TwoDimensional = 1;
    imageScan.slices = 1;
    imageScan.magimage = zeros(128, 128, 1, 'single');
    imageScan.fov = 300;
    imageScan.fov_phase = 300;
    imageScan.fov_3d = 10;
    imageScan.param.OFF_CENTER_FIELD_OF_VIEW_1D = -0.00586;
    imageScan.param.OFF_CENTER_FIELD_OF_VIEW_2D = -0.07031;
    imageScan.param.OFF_CENTER_FIELD_OF_VIEW_3D = 0.01953;
    imageScan.param.OFF_CENTER_FIELD_OF_VIEW_X = -0.00586;
    imageScan.param.OFF_CENTER_FIELD_OF_VIEW_Y = 0.07031;
    imageScan.param.OFF_CENTER_FIELD_OF_VIEW_Z = -0.01953;

    svs = syntheticSVSScan();
    svs.param.OFF_CENTER_FIELD_OF_VIEW_1D = -0.02149;
    svs.param.OFF_CENTER_FIELD_OF_VIEW_2D = -0.02344;
    svs.param.OFF_CENTER_FIELD_OF_VIEW_3D = 0.02734;
    svs.param.OFF_CENTER_FIELD_OF_VIEW_X = -0.02149;
    svs.param.OFF_CENTER_FIELD_OF_VIEW_Y = 0.02344;
    svs.param.OFF_CENTER_FIELD_OF_VIEW_Z = -0.02734;

    imageGeometry = ffc_scan_spatial_geometry(imageScan);
    svsGeometry = ffc_scan_spatial_geometry(svs);
    view = ffc_anatomical_plane_geometry(imageScan, 'Axial', 1, eye(3));
    relativeCentre = svsGeometry.centerPatientMm - ...
        imageGeometry.centerPatientMm;
    displayCentre = [dot(view.basisXPatient, relativeCentre), ...
        dot(view.basisYPatient, relativeCentre)];

    % The scanner prescription is leftward and posterior relative to the
    % spin-echo FOV centre: rightward and downward on a radiological screen.
    verifyEqual(testCase, displayCentre, [15.63 46.87], ...
        'AbsTol', 1e-10);
end


function testVoxelIntersectionProducesClosedOverlay(testCase)
    scan = syntheticImagingScan();
    svs = syntheticSVSScan();
    svsGeometry = ffc_scan_spatial_geometry(svs);
    [sliceIndex, view] = ffc_nearest_anatomical_slice( ...
        scan, 'Axial', svsGeometry.centerPatientMm, eye(3));
    [x, y, details] = ffc_voxel_plane_intersection(svs, view);

    verifyEqual(testCase, sliceIndex, 5);
    verifyTrue(testCase, details.intersects);
    verifyGreaterThanOrEqual(testCase, numel(x), 4);
    verifyEqual(testCase, x(1), x(end), 'AbsTol', 1e-10);
    verifyEqual(testCase, y(1), y(end), 'AbsTol', 1e-10);
    verifyLessThanOrEqual(testCase, max(abs(x)), 10 + 1e-10);
    verifyLessThanOrEqual(testCase, max(abs(y)), 10 + 1e-10);
end


function scan = syntheticImagingScan()
    scan = struct();
    scan.orientation = [1 0 0 0 -1 0];
    scan.fov = 40;
    scan.fov_phase = 60;
    scan.fov_3d = 80;
    scan.thk = 10;
    scan.samples = 4;
    scan.views = 6;
    scan.slices = 8;
    scan.magimage = zeros(4, 6, 8, 'single');
    scan.TwoDimensional = 0;
    scan.is_svs = false;
    scan.param = struct('OFF_CENTER_FIELD_OF_VIEW_1D', 0, ...
        'OFF_CENTER_FIELD_OF_VIEW_2D', 0, ...
        'OFF_CENTER_FIELD_OF_VIEW_3D', 0);
end


function scan = syntheticSVSScan()
    scan = struct();
    scan.orientation = [1 0 0 0 -1 0];
    scan.fov = 20;
    scan.fov_phase = 20;
    scan.fov_3d = 20;
    scan.thk = 20;
    scan.samples = 4;
    scan.views = 1;
    scan.slices = 1;
    scan.magimage = zeros(1, 1, 1, 'single');
    scan.TwoDimensional = 1;
    scan.is_svs = true;
    scan.voxel_size_mm = [20 20 20];
    scan.param = struct('OFF_CENTER_FIELD_OF_VIEW_1D', 0, ...
        'OFF_CENTER_FIELD_OF_VIEW_2D', 0, ...
        'OFF_CENTER_FIELD_OF_VIEW_3D', 0);
end
