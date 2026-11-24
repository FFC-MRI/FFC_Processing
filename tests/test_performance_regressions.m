function tests = test_performance_regressions
% Numerical and state checks for the performance edition. Requires MATLAB.
    tests = functiontests(localfunctions);
end
function setupOnce(testCase)
    root = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.oldPath = path;
    addpath(root); ffc_setup();
end
function teardownOnce(testCase)
    path(testCase.TestData.oldPath);
end

function testBlockedWhiteningMatchesDirectSolve(testCase)
    rng(91);
    K = complex(randn(9,7,3,2,3,4,'single'),randn(9,7,3,2,3,4,'single'));
    info = ffc_noise_calibration(K);
    X = reshape(K,[],4).';
    expected = info.L \ X;
    actual = ffc_apply_whitening(K,info,[4 2],4,17);
    expected = reshape(expected([4 2],:).',[9 7 3 2 3 2]);
    verifyEqual(testCase,actual,expected,'AbsTol',single(3e-6));
end

function testWhiteningEigenFallbackMatchesDirectTransform(testCase)
    K = reshape(complex(single(1:360),single(360:-1:1)),[5 4 3 1 1 6]);
    W = complex(single(eye(6)),single(tril(ones(6),-1)));
    info = struct('psi',eye(6),'p',1,'Winv',W);
    actual = ffc_apply_whitening(K,info,[6 1],6,13);
    X = W * reshape(K,[],6).';
    expected = reshape(X([6 1],:).',[5 4 3 1 1 2]);
    verifyEqual(testCase,actual,expected);
end

function testSingleAllocationPaddingMatchesLegacyConvention(testCase)
    K = reshape(complex(single(1:630)),[5 7 3 2 1 3]);
    target = [8 12 6];
    delta = target - [5 7 3];
    expected = padarray(padarray(K,floor(delta/2),0,'pre'),ceil(delta/2),0,'post');
    verifyEqual(testCase,ffc_center_pad(K,target),expected);
    verifyEqual(testCase,ffc_center_pad(K,[5 7 3]),K);
end

function testBlockedReconstructionMatchesBulk2DAnd3D(testCase)
    for is2D = [true false]
        for whitening = [false true]
            for mode = ["RSSQ" "Adaptive"]
                scan = makeScan(is2D);
                scan.noise_whitening = whitening;
                scan.combination = char(mode);
                scan.fft_size = 14; scan.fft_size_3d = 5;
                scan.multichannel_recon = [4 2];
                scan = buildimages(scan);
                [magnitude,combined,coils] = bulkReference(scan);
                verifyEqual(testCase,scan.magimage,magnitude,'AbsTol',single(5e-5));
                verifyEqual(testCase,scan.complexcombined,combined,'AbsTol',single(5e-5));
                if mode == "RSSQ"
                    verifyEqual(testCase,scan.compleximage,coils,'AbsTol',single(5e-5));
                end
            end
        end
    end
end

function testAcsMapsAverageEveryFieldAndReceiver(testCase)
    rng(5);
    K = complex(randn(9,7,3,2,3,2,'single'),randn(9,7,3,2,3,2,'single'));
    K(:,:,:,:,3,2) = 5*K(:,:,:,:,3,2);
    obj = struct('TwoDimensional',1);
    actual = estimate_maps_from_acs_walsh(K,1:2,obj);
    reference = referenceWalsh(K,true);
    verifyEqual(testCase,actual,reference,'AbsTol',single(2e-5));
end

function testAcsMapsUseSpatial3dTransform(testCase)
    rng(6);
    K = complex(randn(9,7,3,2,2,2,'single'),randn(9,7,3,2,2,2,'single'));
    actual = estimate_maps_from_acs_walsh(K,1:2,struct('TwoDimensional',0));
    expected = referenceWalsh(K,false);
    verifyEqual(testCase,actual,expected,'AbsTol',single(2e-5));
    wrong2D = referenceWalsh(K,true);
    verifyGreaterThan(testCase,norm(actual(:)-wrong2D(:)),0.1);
end

function testBlockedAcsReconstructionMatchesBulk(testCase)
    for is2D = [true false]
        scan = makeScan(is2D); scan.noise_whitening = true;
        scan.combination = 'ACS/Walsh'; scan.fft_size = 14; scan.fft_size_3d = 5;
        scan = buildimages(scan);
        [magnitude,combined] = bulkReference(scan);
        verifyEqual(testCase,scan.magimage,magnitude,'AbsTol',single(5e-5));
        verifyEqual(testCase,scan.complexcombined,combined,'AbsTol',single(5e-5));
    end
end

function testCachedReconstructionPreservesPostprocessing(testCase)
    scan = buildimages(makeScan(true));
    scan.magimage = scan.magimage * 0.75;
    scan.T1Maps = ones(9,7,3); scan.R1Maps = scan.T1Maps;
    scan.param.SimpleIRMapFit = struct('Rsquare',ones(9,7,3));
    scan.param.SimpleIRFit = struct('test',true);
    scan.geomT = [0 1 0;-1 0 0;0 0 1]; scan.rot_theta = 90;
    retained = buildimages(scan,false);
    verifyEqual(testCase,retained.magimage,scan.magimage);
    verifyEqual(testCase,retained.T1Maps,scan.T1Maps);
    verifyEqual(testCase,retained.reconstruction_generation,scan.reconstruction_generation);
    rebuilt = buildimages(scan);
    verifyEmpty(testCase,rebuilt.T1Maps);
    verifyFalse(testCase,isfield(rebuilt.param,'SimpleIRMapFit'));
    verifyEqual(testCase,rebuilt.reconstruction_generation,scan.reconstruction_generation+1);
end

function testSettingsInvalidateReconstructionButReuseCalibration(testCase)
    scan = makeScan(true); scan.noise_whitening = true;
    scan = buildimages(scan);
    oldInfo = scan.whitening_cache;
    generation = scan.reconstruction_generation;
    scan.fft_size = 18;
    verifyFalse(testCase,ffc_reconstruction_valid(scan));
    scan = buildimages(scan,false);
    verifyEqual(testCase,scan.whitening_cache,oldInfo);
    verifyEqual(testCase,scan.reconstruction_generation,generation+1);
    verifyTrue(testCase,ffc_reconstruction_valid(scan));
end

function testRepeatedPreprocessingNeverBakesInOffset(testCase)
    scan = makeScan(true);
    scan.param.OFF_CENTER_FIELD_OF_VIEW_1D = 0.011;
    original = scan.originalcomplexkspace;
    scan = buildimages(scan);
    image = scan.magimage;
    scan.reprocess = false;
    scan = preprocessing(scan);
    verifyEqual(testCase,scan.originalcomplexkspace,original);
    scan = buildimages(scan);
    verifyEqual(testCase,scan.magimage,image);
end

function testRssqCoilRetentionCanBeDisabled(testCase)
    scan = buildimages(makeScan(true));
    magnitude = scan.magimage;
    scan.retain_coil_images = false; scan = buildimages(scan,false);
    verifyEqual(testCase,scan.magimage,magnitude);
    verifyEqual(testCase,size(scan.compleximage,6),1);
end

function testLlrCoverageAndBordersMatchLegacy(testCase)
    rng(12); X = abs(randn(11,13,7,'single'));
    opts = struct('patch',5,'step',3,'niter',2,'tau',0.5,'sigma',0.1);
    actual = llr_denoise_stack(X,opts);
    expected = legacyLlr(X,opts);
    verifyEqual(testCase,actual,expected,'AbsTol',single(3e-6));
    verifyEqual(testCase,actual(:,end,:),X(:,end,:));
end

function testPerSliceGroupingMatchesExactPreviewScope(testCase)
    rng(13); X = abs(randn(11,13,3,2,3,'single'));
    params = struct('patch',5,'step',3,'niter',1,'tau',0.4,'sigma',0.1, ...
        'use_parallel',false,'stack_grouping','per-slice');
    allImages = ffc_mri_filter(X,'LLR',params);
    selectedSlice = ffc_mri_filter(X(:,:,2,:,:),'LLR',params);
    verifyEqual(testCase,allImages(:,:,2,:,:),selectedSlice);
    params = rmfield(params,'stack_grouping');
    legacy = ffc_mri_filter(X,'LLR',params);
    params.stack_grouping = 'all';
    verifyEqual(testCase,ffc_mri_filter(X,'LLR',params),legacy);
end

function testFilterApplyIsNoncumulativeAndClearsMaps(testCase)
    scan = buildimages(makeScan(true));
    scan.T1Maps = ones(9,7,3);
    params = struct('patch',3,'step',2,'niter',1,'tau',0.2,'sigma',0.1,'stack_grouping','per-slice');
    task = struct('filter','LLR','params',params);
    once = ffc_compute_task('filter',scan,task);
    twice = ffc_compute_task('filter',once,task);
    verifyEqual(testCase,twice.magimage,once.magimage);
    verifyEmpty(testCase,once.T1Maps);
    reset = ffc_compute_task('filter',twice,struct('filter','None','params',struct()));
    verifyEqual(testCase,reset.magimage,scan.magimage);
end

function testScalarIrOutputPreservesFullFit(testCase)
    t = [600;10;40;100;250;1200];
    y = -1.8*exp(-t/350)+0.94;
    full = ffc_simple_ir_fit(t,y);
    [T1,R1,r2] = ffc_simple_ir_fit(t,y);
    verifyEqual(testCase,[T1 R1 r2],[full.T1ms full.R1 full.rsquare]);
    verifyEqual(testCase,T1,350,'AbsTol',1e-3);
    y(3) = NaN;
    full = ffc_simple_ir_fit(t,y);
    [T1,R1,r2] = ffc_simple_ir_fit(t,y);
    verifyEqual(testCase,[T1 R1 r2],[full.T1ms full.R1 full.rsquare]);
end

function testMapKernelFitsOnlySelectedVoxels(testCase)
    scan = makeScan(true);
    t = [10 40 100 250 600 1200];
    scan.timepoints = t; scan.fieldpoints = 200;
    scan.magimage = repmat(reshape(single(abs(-1.8*exp(-t/350)+0.94)),1,1,1,6),[3 4 2 1]);
    scan.mask = false(3,4,2); scan.mask(2,3,1) = true;
    scan = ffc_simple_t1map(scan,struct('use_parallel',false));
    verifyEqual(testCase,nnz(isfinite(scan.T1Maps)),1);
    verifyEqual(testCase,scan.T1Maps(2,3,1),350,'AbsTol',0.1);
end

function scan = makeScan(is2D)
    rng(77);
    shape = [9 7 3]; nTimes = 2; nFields = 3; nCoils = 4;
    raw = complex(randn([shape nTimes*nFields nCoils],'single'),randn([shape nTimes*nFields nCoils],'single'));
    c = struct('NUMBER_OF_AVERAGES',1,'ACCU_DIM',1,'FIELD_OF_VIEW',0.09, ...
        'FIELD_OF_VIEW_PHASE',0.07,'FIELD_OF_VIEW_3D',0.03,'SLICE_THICKNESS',0.01, ...
        'IMAGE_ORIENTATION_SUBJECT',[1 0 0 0 -1 0],'MULTI_PLANAR_EXCITATION',is2D, ...
        'USER_MATRIX_DIMENSION_2D',shape(2),'ACQUISITION_MATRIX_DIMENSION_2D',shape(2), ...
        'TRANSFORM_PLUGIN','Sequential4DBackAndForth','FCI_NUMBER_OF_TEVO',nTimes, ...
        'FCI_LIST_OF_T1_ESTIMATES',0.25*ones(1,nFields),'FCI_T1_FACTOR_LOWER',0.05,'FCI_T1_FACTOR_UPPER',8, ...
        'FCI_T1_LOG_DIST',true);
    c.FCI_LIST_OF_BEVO = 200./(1:nFields);
    par = struct('cameleon',c,'md1d',shape(1),'md2d',shape(2),'md3d',shape(3), ...
        'md4d',nTimes*nFields,'ncoils',nCoils,'navg',1,'te',0.01,'fa',90, ...
        'serieName','GRADIENT_ECHO','seqName','GRADIENT_ECHO');
    scan = ImageReconCore(struct('data',raw,'par',par),'.dat',1);
    scan.reprocess = false;
    scan = preprocessing(scan);
    scan.combination = 'RSSQ'; scan.biascorrect = 'Off'; scan.noise_whitening = false;
    scan.fft_size = shape(1); scan.fft_size_3d = shape(3);
    scan.window_function = 'None'; scan.window_size = 1;
    scan.denoise_filter = 'None'; scan.denoise_params = struct();
end

function [magnitude,combined,coils] = bulkReference(scan)
    K = scan.complexkspace; nc = size(K,6); info = scan.whitening_info;
    if isfield(info,'psi') && ~isempty(info.psi)
        X = reshape(K,[],nc).'; X = info.L \ X;
        K = reshape(X.',size(K));
    end
    K = K(:,:,:,:,:,scan.selected_coils_original);
    target = [scan.fft_size round(size(K,2)*scan.fft_size/size(K,1)) size(K,3)];
    if scan.TwoDimensional ~= 1, target(3) = scan.fft_size_3d; end
    delta = target - [size(K,1) size(K,2) size(K,3)];
    K = windowkspace(K,scan.window_size,scan.window_function);
    K = padarray(padarray(K,floor(delta/2),0,'pre'),ceil(delta/2),0,'post');
    if scan.TwoDimensional == 1, coils = ifft2c(K); else, coils = ifft3c(K); end
    if strcmp(scan.combination,'Adaptive')
        ref = 'signalcoil'; if scan.noise_whitening, ref = 'rss'; end
        combined = combine_adaptive_smoothed(coils,struct('coils',1:size(K,6), ...
            'sigma',4,'maskFrac',0.1,'referenceMode',ref));
        magnitude = abs(combined);
    elseif strcmp(scan.combination,'ACS/Walsh')
        maps = referenceWalsh(K,scan.TwoDimensional == 1);
        combined = combine_with_maps(coils,maps,1:size(K,6),scan);
        magnitude = abs(combined);
    else
        magnitude = rssq(coils,6); combined = sum(coils,6);
    end
end

function maps = referenceWalsh(K,is2D)
    [X,Y,S,~,~,C] = size(K);
    cx = floor(X/2)+1; cy = floor(Y/2)+1;
    xr = max(1,cx-12):min(X,cx+11); yr = max(1,cy-12):min(Y,cy+11);
    avg = zeros(X,Y,S,1,1,C,'like',K);
    avg(xr,yr,:,1,1,:) = mean(mean(K(xr,yr,:,:,:,:),5,'omitnan'),4,'omitnan');
    if is2D, low = ifft2c(avg); else, low = ifft3c(avg); end
    maps = zeros(size(avg),'like',K);
    for s = 1:S
        L = reshape(low(:,:,s,1,1,:),X,Y,C);
        rss = sqrt(sum(abs(L).^2,3)); mask = rss > max(eps(class(rss)),0.02*prctile(rss(:),95));
        mf = imgaussfilt(single(mask),3,'FilterSize',19);
        den = imgaussfilt(mf,3,'FilterSize',19)+eps('single');
        R = zeros(X,Y,C,C,'like',L);
        for i = 1:C
            for j = 1:C
                v = L(:,:,i).*conj(L(:,:,j));
                R(:,:,i,j) = complex(imgaussfilt(real(v).*mf,3,'FilterSize',19), ...
                    imgaussfilt(imag(v).*mf,3,'FilterSize',19))./den;
            end
        end
        smap = zeros(X,Y,C,'like',L);
        for x = 1:X
            for y = 1:Y
                if ~mask(x,y), continue; end
                M = squeeze(R(x,y,:,:)); [V,D] = eig((M+M')/2,'vector');
                [~,ix] = max(real(D)); v = V(:,ix); [~,ref] = max(abs(v));
                smap(x,y,:) = v*exp(-1i*angle(v(ref)));
            end
        end
        norm = sqrt(sum(abs(smap).^2,3)); smap = smap./(norm+eps(class(norm)));
        maps(:,:,s,1,1,:) = reshape(smap,[X Y 1 1 1 C]);
    end
end

function Y = legacyLlr(X,opts)
    [H,W,N] = size(X); Y = X; p = opts.patch;
    for pass = 1:opts.niter
        acc = zeros(H,W,N,'single'); wgt = acc;
        for i = 1:opts.step:H-p+1
            for j = 1:opts.step:W-p+1
                patch = Y(i:i+p-1,j:j+p-1,:); P = reshape(patch,p*p,N);
                P = P-mean(P,1,'omitnan'); [U,S,V] = svd(P,'econ');
                sing = max(diag(S)-opts.tau*opts.sigma*sqrt(max(p*p,N)),0);
                P = U*(diag(sing)*V')+mean(reshape(patch,p*p,N),1,'omitnan');
                acc(i:i+p-1,j:j+p-1,:) = acc(i:i+p-1,j:j+p-1,:)+reshape(P,p,p,N);
                wgt(i:i+p-1,j:j+p-1,:) = wgt(i:i+p-1,j:j+p-1,:)+1;
            end
        end
        idx = wgt>0; Y(idx) = acc(idx)./wgt(idx); Y = max(Y,0);
    end
end

function testViewerBatchReadyAndSwitchingDoesNotRebuild(testCase)
    input = {makeScan(true),makeScan(false)};
    [scans,names,failures] = ffc_preprocess_scans(input, ...
        struct('reconstruct',true,'reprocess',false));
    verifyEmpty(testCase,failures); verifyNumElements(testCase,names,2);
    generation = cellfun(@(x) x.reconstruction_generation,scans);
    for index = [1 2 1 2 1]
        scans{index} = ffc_normalize_viewer_scan(scans{index});
        verifyTrue(testCase,ffc_reconstruction_valid(scans{index}));
        scans{index} = buildimages(scans{index},false);
    end
    verifyEqual(testCase,cellfun(@(x) x.reconstruction_generation,scans),generation);
    verifyFalse(testCase,scans{1}.retain_coil_images);
end

function testAnalysisMetadataDoesNotInvalidateOrEraseResults(testCase)
    scan = buildimages(ffc_normalize_viewer_scan(makeScan(true)),false);
    generation = scan.reconstruction_generation;
    scan.param.CVIDispersionMask = struct('mask',true(9,7),'registration',eye(4));
    scan.param.T1DispersionModelFit = struct('model','power law','R1',1:3);
    scan.param.LastEnhancedDICOMExport = struct('folder','example');
    scan.T1FitResults = struct('T1ms',[100 200 300]);
    scan.mask = true(9,7); scan.T1Maps = ones(9,7,3);
    scan.registration_tform = struct('test',true);
    before = scan;
    scan = buildimages(scan,false);
    verifyEqual(testCase,scan.reconstruction_generation,generation);
    verifyEqual(testCase,scan.T1FitResults,before.T1FitResults);
    verifyEqual(testCase,scan.mask,before.mask);
    verifyEqual(testCase,scan.T1Maps,before.T1Maps);
    verifyEqual(testCase,scan.registration_tform,before.registration_tform);
    rebuilt = buildimages(scan,true);
    verifyEmpty(testCase,rebuilt.T1FitResults);
    verifyFalse(testCase,isfield(rebuilt.param,'T1DispersionModelFit'));
    verifyFalse(testCase,isfield(rebuilt.param,'CVIDispersionMask'));
end

function testRealSettingsAndOffsetsStillInvalidate(testCase)
    scan = buildimages(ffc_normalize_viewer_scan(makeScan(true)),false);
    settings = {'fft_size','window_size','noise_whitening','multichannel_recon'};
    values = {18,0.8,true,[2 1]};
    for k = 1:numel(settings)
        changed = scan; changed.(settings{k}) = values{k};
        verifyFalse(testCase,ffc_reconstruction_valid(changed));
    end
    for name = {'OFF_CENTER_FIELD_OF_VIEW_1D','OFF_CENTER_FIELD_OF_VIEW_2D', ...
            'OFF_CENTER_FIELD_OF_VIEW_3D','FIELD_OF_VIEW','FIELD_OF_VIEW_PHASE','FIELD_OF_VIEW_3D'}
        changed = scan; changed.param.(name{1}) = 0.013;
        verifyFalse(testCase,ffc_reconstruction_valid(changed));
    end
    changed = scan; changed.param.OFF_CENTER_FIELD_OF_VIEW_EFF = [0.01 0 0];
    verifyFalse(testCase,ffc_reconstruction_valid(changed));
end

function testViewerDefaultsAreIdempotent(testCase)
    scan = makeScan(true);
    scan.multichannel_recon = logical([1 0 1 0]);
    scan.window_function = 'none'; scan.fft_size_3d = [];
    scan = ffc_normalize_viewer_scan(scan);
    key = ffc_reconstruction_key(scan);
    scan = ffc_normalize_viewer_scan(scan);
    verifyEqual(testCase,ffc_reconstruction_key(scan),key);
    verifyEqual(testCase,scan.multichannel_recon,[1 3]);
    verifyEqual(testCase,scan.window_function,'None');
end

function testDirectSlabMatchesVolumeExtraction(testCase)
    for is2D = [true false]
        scan = buildimages(makeScan(is2D));
        for name = {'Axial','Coronal','Sagittal'}
            for transform = {eye(3),[0 -1 0;1 0 0;0 0 1],[-1 0 0;0 1 0;0 0 1]}
                volume = scan.magimage(:,:,:,2,3);
                expected = ffc_extract_anatomical_plane(volume,scan,name{1},2,transform{1});
                actual = ffc_extract_anatomical_plane(scan.magimage,scan,name{1},2,transform{1},2,3);
                verifyEqual(testCase,actual,expected);
            end
        end
    end
end

function testDeferredPhaseMatchesStoredPhase(testCase)
    full = buildimages(makeScan(false));
    lean = buildimages(ffc_normalize_viewer_scan(makeScan(false)),false);
    verifyEmpty(testCase,lean.phaseimage);
    for name = {'Axial','Coronal','Sagittal'}
        expected = ffc_extract_anatomical_plane(full.phaseimage,full,name{1},2,eye(3),2,3);
        actual = angle(ffc_extract_anatomical_plane(lean.complexcombined,lean,name{1},2,eye(3),2,3));
        verifyEqual(testCase,actual,expected);
    end
end

function testSensitivityCacheAndInvalidation(testCase)
    scan = makeScan(false); scan.combination = 'ACS/Walsh';
    scan = buildimages(scan,false);
    key = scan.sensitivity_cache_key; maps = scan.sensitivity_cache;
    scan.biascorrect = 'ACS';
    reused = buildimages(scan,false);
    forced = buildimages(scan,true);
    verifyEqual(testCase,reused.sensitivity_cache_key,key);
    verifyEqual(testCase,reused.sensitivity_cache,maps);
    verifyEqual(testCase,reused.magimage,forced.magimage,'AbsTol',single(5e-5));
    reused.multichannel_recon = [2 1]; reused = buildimages(reused,false);
    verifyNotEqual(testCase,reused.sensitivity_cache_key,key);
    verifyEqual(testCase,size(reused.sensitivity_cache,6),2);
end

function testBatchCancellationAndFailureIsolation(testCase)
    good = makeScan(true); bad = good;
    bad.fft_size = 1;
    [scans,names,failures] = ffc_preprocess_scans({good,bad,good},struct('reconstruct',true));
    verifyNumElements(testCase,scans,2); verifyNumElements(testCase,names,2);
    verifyEqual(testCase,failures.sequenceIndex,2);
    verifyTrue(testCase,all(cellfun(@ffc_reconstruction_valid,scans)));
    verifyError(testCase,@() ffc_preprocess_scans({good,good}, ...
        struct('reconstruct',true,'progress',@(i,n) i < 2)), 'ffc:LoadCancelled');
end

function testDiskCacheRoundTripAndSourceChanges(testCase)
    root = tempname; mkdir(root); cache = fullfile(root,'cache');
    names = {'DiskCacheDirectory','DiskCacheEnabled'};
    old = cell(1,2); existed = false(1,2);
    for k = 1:2
        existed(k) = ispref('FFCProcessing',names{k});
        if existed(k), old{k} = getpref('FFCProcessing',names{k}); end
    end
    cleanup = onCleanup(@restore); %#ok<NASGU>
    setpref('FFCProcessing','DiskCacheDirectory',cache);
    setpref('FFCProcessing','DiskCacheEnabled',true);
    writeBytes(fullfile(root,'data.dat'),uint8(1:20));
    writeBytes(fullfile(root,'header.xml'),uint8('<header>one</header>'));
    writeBytes(fullfile(root,'Serie.xml'),uint8('<serie>one</serie>'));
    scan = buildimages(ffc_normalize_viewer_scan(makeScan(true)),false);
    scan.disk_cache_source = root;
    scan.disk_cache_identity = ffc_disk_cache('identity',root);
    verifyTrue(testCase,ffc_disk_cache('save',root,scan));
    restored = ffc_disk_cache('load',root);
    verifyEqual(testCase,restored.magimage,scan.magimage);
    verifyEqual(testCase,restored.reconstruction_generation,scan.reconstruction_generation);
    % Metadata changes invalidate even when the byte length stays the same.
    writeBytes(fullfile(root,'Serie.xml'),uint8('<serie>two</serie>'));
    verifyEmpty(testCase,ffc_disk_cache('load',root));
    writeBytes(fullfile(root,'Serie.xml'),uint8('<serie>one</serie>'));
    writeBytes(fullfile(root,'header.xml'),uint8('<header>two</header>'));
    verifyEmpty(testCase,ffc_disk_cache('load',root));
    writeBytes(fullfile(root,'header.xml'),uint8('<header>one</header>'));
    writeBytes(fullfile(root,'data.dat'),uint8(1:21));
    verifyEmpty(testCase,ffc_disk_cache('load',root));
    verifyEmpty(testCase,ffc_disk_cache('save',root,scan));
    function restore()
        for j = 1:2
            if existed(j), setpref('FFCProcessing',names{j},old{j});
            elseif ispref('FFCProcessing',names{j}), rmpref('FFCProcessing',names{j}); end
        end
        if isfolder(root), rmdir(root,'s'); end
    end
end
function writeBytes(file,bytes)
    fid = fopen(file,'wb'); cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid,bytes,'uint8');
end
