function tests = test_network_pipeline
    tests = functiontests(localfunctions);
end
function setupOnce(testCase)
    testCase.TestData.oldPath = path;
    addpath(fileparts(fileparts(mfilename('fullpath')))); ffc_setup();
end
function teardownOnce(testCase)
    path(testCase.TestData.oldPath);
end
function setup(testCase)
    names = {'ParallelLoading','LoadWorkers','DiskCacheEnabled'};
    old = cell(size(names)); present = false(size(names));
    for k=1:numel(names)
        present(k)=ispref('FFCProcessing',names{k});
        if present(k), old{k}=getpref('FFCProcessing',names{k}); end
    end
    testCase.TestData.names=names; testCase.TestData.old=old; testCase.TestData.present=present;
    setpref('FFCProcessing','DiskCacheEnabled',false);
    setpref('FFCProcessing','LoadWorkers',2);
    setpref('FFCProcessing','ParallelLoading',false);
end
function teardown(testCase)
    for k=1:numel(testCase.TestData.names)
        name=testCase.TestData.names{k};
        if testCase.TestData.present(k), setpref('FFCProcessing',name,testCase.TestData.old{k});
        elseif ispref('FFCProcessing',name), rmpref('FFCProcessing',name); end
    end
end
function testSerialFallbackOrderAndFailureIsolation(testCase)
    [scans,report,failures]=ffc_load_raw_pipeline({'1','bad','3'},'',struct('reader',@reader));
    verifyEqual(testCase,numel(scans),2);
    verifyEqual(testCase,failures.sequenceIndex,2);
    verifyEqual(testCase,scans{1}.name,'1'); verifyEqual(testCase,scans{2}.name,'3');
    verifyTrue(testCase,all(cellfun(@ffc_reconstruction_valid,scans)));
    verifyEqual(testCase,report.workers,0);
    verifyEqual(testCase,report.successCount,2);
end
function testCancellationBeforeReading(testCase)
    verifyError(testCase,@() ffc_load_raw_pipeline({'1','2'},'', ...
        struct('reader',@reader,'loadProgress',@(state) false)), 'ffc:LoadCancelled');
end
function testWarmCacheBypassesWorkers(testCase)
    setpref('FFCProcessing','ParallelLoading',true);
    [scans,report,failures]=ffc_load_raw_pipeline({'1','2'},'',struct('reader',@warmReader));
    verifyEmpty(testCase,failures); verifyEqual(testCase,report.workers,0);
    verifyEqual(testCase,report.maxInFlight,0);
    verifyTrue(testCase,all([report.sequence.cacheHit]));
    verifyEqual(testCase,scans{1}.reconstruction_generation,uint64(1));
end
function testParallelMatchesSerialAndQueueBound(testCase)
    assumeTrue(testCase,license('test','Distrib_Computing_Toolbox'));
    pool=gcp('nocreate');
    assumeFalse(testCase,isa(pool,'parallel.ThreadPool'));
    expected=ffc_load_raw_pipeline({'1','2','3','4'},'',struct('reader',@reader));
    setpref('FFCProcessing','ParallelLoading',true);
    [actual,report,failures]=ffc_load_raw_pipeline({'1','2','3','4'},'', ...
        struct('reader',@reader,'loadProgress',@checkState));
    verifyEmpty(testCase,failures);
    verifyGreaterThan(testCase,report.workers,0);
    verifyLessThanOrEqual(testCase,report.maxInFlight,2);
    for k=1:4
        verifyEqual(testCase,actual{k}.name,expected{k}.name);
        verifyEqual(testCase,actual{k}.magimage,expected{k}.magimage,'AbsTol',single(5e-5));
    end
    function proceed=checkState(state)
        verifyLessThanOrEqual(testCase,state.active,2);
        proceed=true;
    end
end
function testCancelWithJobInFlight(testCase)
    assumeTrue(testCase,license('test','Distrib_Computing_Toolbox'));
    pool=gcp('nocreate'); assumeFalse(testCase,isa(pool,'parallel.ThreadPool'));
    setpref('FFCProcessing','ParallelLoading',true);
    verifyError(testCase,@() ffc_load_raw_pipeline({'1','2','3'},'', ...
        struct('reader',@reader,'loadProgress',@(state) state.active==0)), 'ffc:LoadCancelled');
end
function [scan,timing]=reader(folder,~,~)
    if strcmp(folder,'bad'), error('fixture:ReadFailed','Expected reader failure'); end
    scan=makeScan(true); scan.name=folder;
    scan.rawdata=scan.rawdata*str2double(folder);
    scan.originalcomplexkspace=[]; scan.disk_cache_source='';
    timing=struct('readSeconds',0,'cacheLookupSeconds',0,'cacheHit',false);
end
function [scan,timing]=warmReader(folder,root,reprocess)
    [scan,timing]=reader(folder,root,reprocess);
    result=ffc_prepare_loaded_scan(scan,struct()); scan=result.scan;
    timing.cacheHit=true;
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
