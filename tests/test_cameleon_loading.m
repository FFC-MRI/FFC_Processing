function tests = test_cameleon_loading
% Regression tests for stored averages, receiver dimensions and failed loads.
    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.originalPath = path;
    addpath(projectRoot);
    addpath(fullfile(projectRoot, 'misc_functions'));
    addpath(fullfile(projectRoot, 'kspace_handling'));
    addpath(fullfile(projectRoot, 'image_handling'));
end


function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end


function testAcquisitionAveragesDoNotMultiplyStoredData(testCase)
    for is2D = [true false]
        for nReceivers = [1 13]
            for nAverages = [1 6]
                scan = syntheticScan([8 6 3], nReceivers, nAverages, is2D);
                expected = reshape(scan.rawdata, [8 6 3 1 1 nReceivers]);

                verifyEqual(testCase, scan.averages, 1);
                verifyEqual(testCase, scan.param.NUMBER_OF_AVERAGES, nAverages);
                scan = preprocessing(scan);

                verifyEqual(testCase, size(scan.originalcomplexkspace, 1:6), ...
                    [8 6 3 1 1 nReceivers]);
                verifyEqual(testCase, scan.originalcomplexkspace, expected);
                verifyEqual(testCase, scan.n_receivers, nReceivers);
                verifyEqual(testCase, scan.complexkspace, expected);
                scan.reprocess = true;
                scan = preprocessing(scan);
                verifyEqual(testCase, scan.originalcomplexkspace, expected);
            end
        end
    end
end


function testReported128By128By15ThirteenReceiverScan(testCase)
    % The supplied GRADIENT_ECHO scan stores 3,194,880 complex samples,
    % despite NUMBER_OF_AVERAGES=6. Its byte count is 25,559,040 (float32 IQ).
    scan = syntheticScan([128 128 15], 13, 6, true);
    raw = scan.rawdata;
    verifyEqual(testCase, 8 * numel(raw), 25559040);

    scan = preprocessing(scan);

    verifyEqual(testCase, size(scan.originalcomplexkspace, 1:6), ...
        [128 128 15 1 1 13]);
    verifyEqual(testCase, scan.originalcomplexkspace, ...
        reshape(raw, [128 128 15 1 1 13]));
end


function testImagingReconstructionKeepsAllThirteenReceivers(testCase)
    for is2D = [true false]
        scan = syntheticScan([8 6 3], 13, 6, is2D);
        scan = preprocessing(scan);
        original = scan.originalcomplexkspace;
        scan.multichannel_recon = 1:13;
        scan.noise_whitening = false;
        scan.combination = 'RSSQ';
        scan.biascorrect = 'Off';

        scan = buildimages(scan);
        if is2D
            expected = ifft2c(original);
        else
            expected = ifft3c(original);
        end

        verifyEqual(testCase, scan.selected_coils_original, 1:13);
        verifyEqual(testCase, size(scan.compleximage, 6), 13);
        verifyEqual(testCase, scan.originalcomplexkspace, original);
        verifyEqual(testCase, scan.compleximage, expected, 'AbsTol', single(1e-5));
        verifyEqual(testCase, scan.magimage, rssq(expected, 6), 'AbsTol', single(1e-5));
        verifyTrue(testCase, all(isfinite(scan.magimage(:))));
    end
end


function testFieldAndTimeOrderIsPreserved(testCase)
    scan = syntheticScan([8 6 3], 13, 6, true, false, 2, 3);
    expected = reshape(scan.rawdata, [8 6 3 2 3 13]);

    scan = preprocessing(scan);

    verifyEqual(testCase, [scan.n_timepoints scan.n_fieldpoints], [2 3]);
    verifyEqual(testCase, scan.originalcomplexkspace, expected);
end


function testSVSStillUsesMagnitudeBeforeMeanWithSixAccumulations(testCase)
    scan = syntheticScan([8 1 1], 1, 6, true, true, 2, 3);
    raw = scan.rawdata;
    raw(2:2:end, :, :, :, :) = -raw(2:2:end, :, :, :, :);
    scan.rawdata = raw;
    scan = preprocessing(scan);
    scan = buildimages(scan);

    verifyEqual(testCase, scan.averages, 1);
    verifyEqual(testCase, scan.param.NUMBER_OF_AVERAGES, 6);
    verifyEqual(testCase, scan.magimage, ...
        reshape(mean(abs(raw), 1), [1 1 1 2 3]), 'AbsTol', single(1e-5));
end


function testExplicitStoredAverageDimensionStillGetsAveraged(testCase)
    % Legacy formats really can store separate averages. Do not change the
    % generic preprocessing mean or globally force every format to one.
    scan = syntheticScan([8 6 3], 13, 1, true, false, 2, 3);
    first = reshape(scan.rawdata, [8 6 1 3 2 3 1 13]);
    stored = cat(7, first, first + complex(single(20), single(-10)));
    scan.averages = 2;
    scan.rawdata = reshape(stored, 8, 6, []);
    expected = reshape(mean(stored, 7), [8 6 3 2 3 13]);

    scan = preprocessing(scan);

    verifyEqual(testCase, scan.originalcomplexkspace, expected);
end


function testMismatchedSampleCountHasDimensionDiagnostic(testCase)
    scan = syntheticScan([8 6 3], 13, 6, true);
    scan.samples = 9;
    verifyError(testCase, @() preprocessing(scan), ...
        'ImageReconCore:ImagingDataSizeMismatch');
end


function testInvalidStoredDimensionHasDiagnostic(testCase)
    scan = syntheticScan([8 6 3], 13, 6, true);
    scan.averages = 0;
    verifyError(testCase, @() preprocessing(scan), ...
        'ImageReconCore:InvalidImagingDimensions');
end


function testMissingPreprocessingDoesNotReportOneAvailableCoil(testCase)
    scan = syntheticScan([8 6 3], 13, 6, true);
    scan.multichannel_recon = 1:13;
    verifyError(testCase, @() buildimages(scan), ...
        'ImageReconCore:NoPreprocessedKspace');
end


function testRealOutOfRangeReceiverSelectionStillFails(testCase)
    scan = syntheticScan([8 6 3], 1, 6, true);
    scan = preprocessing(scan);
    scan.multichannel_recon = 1:13;
    verifyError(testCase, @() buildimages(scan), ...
        'buildimages:coilIndexOutOfRange');
end


function testMixedBatchSkipsFailuresAndKeepsNamesAligned(testCase)
    good = syntheticScan([8 6 3], 13, 6, true);
    bad = good;
    bad.samples = 9;
    secondGood = syntheticScan([8 6 3], 1, 1, false);
    secondGood.sequence = 'Second valid scan';
    options = struct('backgroundselect', false, 'reprocess', true, 'checkfit', false);

    [scans, names, failures] = ffc_preprocess_scans( ...
        {bad, good, bad, secondGood}, options);

    verifyEqual(testCase, numel(scans), 2);
    verifyEqual(testCase, names, {good.sequence, secondGood.sequence});
    verifyEqual(testCase, [failures.sequenceIndex], [1 3]);
    verifyEqual(testCase, failures(1).exception.identifier, ...
        'ImageReconCore:ImagingDataSizeMismatch');
    verifyNotEmpty(testCase, failures(1).exception.stack);
    verifyTrue(testCase, scans{1}.reprocess);
    verifyFalse(testCase, scans{1}.backgroundselect);
    verifyFalse(testCase, scans{1}.checkfit);
    verifyEqual(testCase, size(scans{1}.originalcomplexkspace, 6), 13);
    verifyEqual(testCase, size(scans{2}.originalcomplexkspace, 6), 1);
    verifyEmpty(testCase,scans{1}.T1Maps); % maps are now allocated on demand
    verifyEmpty(testCase,scans{1}.R1Maps);
end


function testAllSuccessfulBatchKeepsEveryScan(testCase)
    first = syntheticScan([8 6 3], 13, 6, true);
    second = syntheticScan([8 6 3], 1, 1, false);
    second.sequence = 'Second valid scan';

    [scans, names, failures] = ffc_preprocess_scans({first, second});

    verifyEqual(testCase, numel(scans), 2);
    verifyEqual(testCase, names, {first.sequence, second.sequence});
    verifyEqual(testCase, scans{1}.sequence, first.sequence);
    verifyEqual(testCase, scans{2}.sequence, second.sequence);
    verifyEmpty(testCase, failures);
end


function testAllFailedBatchReturnsNoDisplayableScans(testCase)
    bad = syntheticScan([8 6 3], 13, 6, true);
    bad.samples = 9;

    [scans, names, failures] = ffc_preprocess_scans({bad, bad});

    verifyEmpty(testCase, scans);
    verifyEmpty(testCase, names);
    verifyEqual(testCase, [failures.sequenceIndex], [1 2]);
end


function testEmptyBatchReturnsEmptyResults(testCase)
    [scans, names, failures] = ffc_preprocess_scans({});
    verifyEmpty(testCase, scans);
    verifyEmpty(testCase, names);
    verifyEmpty(testCase, failures);
end


function scan = syntheticScan(matrixSize, nReceivers, nAverages, is2D, isSVS, nTimes, nFields)
    if nargin < 5
        isSVS = false;
    end
    if nargin < 6
        nTimes = 1;
    end
    if nargin < 7
        nFields = 1;
    end
    nExperiments = nTimes * nFields;
    shape = [matrixSize nExperiments nReceivers];
    values = reshape(single(1:prod(shape)), shape);
    raw = complex(values, -values);
    cameleon = struct( ...
        'NUMBER_OF_AVERAGES', nAverages, ...
        'ACCU_DIM', 1, ...
        'FIELD_OF_VIEW', 0.080, ...
        'FIELD_OF_VIEW_PHASE', 0.060, ...
        'FIELD_OF_VIEW_3D', 0.030, ...
        'SLICE_THICKNESS', 0.010, ...
        'IMAGE_ORIENTATION_SUBJECT', [1 0 0 0 -1 0], ...
        'MULTI_PLANAR_EXCITATION', is2D, ...
        'USER_MATRIX_DIMENSION_2D', matrixSize(2), ...
        'ACQUISITION_MATRIX_DIMENSION_2D', matrixSize(2), ...
        'TRANSFORM_PLUGIN', 'Sequential4DBackAndForth');
    if nExperiments > 1
        cameleon.FCI_LIST_OF_BEVO = 200 ./ (1:nFields);
        cameleon.FCI_NUMBER_OF_TEVO = nTimes;
        cameleon.FCI_LIST_OF_T1_ESTIMATES = 0.25 .* ones(1, nFields);
        cameleon.FCI_T1_FACTOR_LOWER = 0.05;
        cameleon.FCI_T1_FACTOR_UPPER = 8;
        cameleon.FCI_T1_LOG_DIST = true;
    end
    sequenceName = 'GRADIENT_ECHO';
    if isSVS
        sequenceName = 'Field-Cycling SVS';
    end
    par = struct('cameleon', cameleon, ...
        'md1d', matrixSize(1), 'md2d', matrixSize(2), ...
        'md3d', matrixSize(3), 'md4d', nExperiments, ...
        'ncoils', nReceivers, 'navg', nAverages, ...
        'te', 0.01, 'fa', 90, 'serieName', sequenceName, 'seqName', sequenceName);
    scan = ImageReconCore(struct('data', raw, 'par', par), '.dat', 1);
    scan.reprocess = false;
end
