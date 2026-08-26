function tests = test_reconstruction_geometry
% Regression tests for centred FFTs, receiver isolation and whitening.

    tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.originalPath = path;
    addpath(fullfile(projectRoot, 'misc_functions'));
    addpath(fullfile(projectRoot, 'kspace_handling'));
end


function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end


function testOddCentredIfftHasNoPhaseRamp(testCase)
    kspace = complex(zeros(5, 7, 'single'));
    kspace(3, 4) = 1;

    image = ifft2c(kspace);

    verifyLessThan(testCase, max(abs(angle(image(:)))), 1e-6);
    verifyLessThan(testCase, max(abs(abs(image(:)) - abs(image(1)))), 1e-6);
end


function testIfft3cDoesNotTransformReceiverDimension(testCase)
    kspace = complex(zeros(7, 5, 3, 1, 1, 2, 'single'));
    kspace(4, 3, 2, 1, 1, 1) = 1;
    kspace(4, 3, 2, 1, 1, 2) = 2i;

    image = ifft3c(kspace);
    receiver1 = image(:,:,:,:,:,1);
    receiver2 = image(:,:,:,:,:,2);

    verifyLessThan(testCase, ...
        max(abs(receiver2(:) - 2i*receiver1(:))), 1e-6);
end


function testWhiteningDecorrelatesSelectedChannels(testCase)
    rng(7);
    matrixSize = [24, 20, 5];
    common = complex(randn(matrixSize, 'single'), randn(matrixSize, 'single'));
    independent1 = complex(randn(matrixSize, 'single'), randn(matrixSize, 'single'));
    independent2 = complex(randn(matrixSize, 'single'), randn(matrixSize, 'single'));

    kspace = complex(zeros([matrixSize, 2], 'single'));
    kspace(:,:,:,1) = common + 0.5*independent1;
    kspace(:,:,:,2) = 0.8*common + independent2;

    whitened = noise_whiten(kspace, struct( ...
        'outer_frac', 0.3, 'max_pages', matrixSize(3), ...
        'max_samples', Inf, 'robust_clip', false));

    samples = reshape(permute(whitened, [4 1 2 3]), 2, []);
    samples = samples - mean(samples, 2);
    covariance = (samples*samples')/(size(samples, 2)-1);

    verifyEqual(testCase, real(diag(covariance)), ones(2, 1), ...
        'AbsTol', 0.08);
    verifyLessThan(testCase, abs(covariance(1, 2)), 0.08);
end
