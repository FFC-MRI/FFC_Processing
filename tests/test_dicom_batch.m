function tests = test_dicom_batch
    tests = functiontests(localfunctions);
end
function setupOnce(tc)
    tc.TestData.oldPath = path;
    addpath(fileparts(fileparts(mfilename('fullpath')))); ffc_setup();
end
function teardownOnce(tc)
    path(tc.TestData.oldPath);
end
function testEveryParticipantAndUniqueOutputs(tc)
    a = scanMeta('A','Same Name','20261110');
    b = scanMeta('B','Same Name','20261110');
    later = scanMeta('A','Same Name','20261111');
    [p,skipped] = ffc_plan_dicom_batch({a,a,b,later},tempdir);
    verifyEmpty(tc,skipped); verifyEqual(tc,[p.index],1:4);
    verifyEqual(tc,numel(unique({p.filename})),4);
    verifyEqual(tc,p(1).options.StudyInstanceUID,p(2).options.StudyInstanceUID);
    verifyNotEqual(tc,p(1).options.StudyInstanceUID,p(3).options.StudyInstanceUID);
    verifyNotEqual(tc,p(1).options.StudyInstanceUID,p(4).options.StudyInstanceUID);
    verifyNotEqual(tc,p(1).options.SeriesInstanceUID,p(2).options.SeriesInstanceUID);
    verifyEqual(tc,p(3).options.ResolvedMetadata.PatientID,'B');
end
function testMissingDemographicsFilledOnlyWithinExplicitIdentity(tc)
    a=scanMeta('A','Alpha','20261110'); b=scanMeta('B','Beta','20261110');
    a.source_metadata.PatientBirthDate='19800101';
    b.source_metadata.PatientBirthDate='19900202';
    partial=a; partial.source_metadata.PatientBirthDate='';
    [p,skipped]=ffc_plan_dicom_batch({a,b,partial},tempdir);
    verifyEmpty(tc,skipped);
    verifyEqual(tc,p(3).options.ResolvedMetadata.PatientBirthDate,'19800101');
    conflicting=a; conflicting.source_metadata.PatientBirthDate='19700101';
    p=ffc_plan_dicom_batch({a,conflicting,partial},tempdir);
    verifyEmpty(tc,p(3).options.ResolvedMetadata.PatientBirthDate);
end
function testUnknownIdentityAndSvsAreReported(tc)
    unknown=scanMeta('','','20261110'); svs=scanMeta('A','Alpha','20261110'); svs.is_svs=true;
    [p,skipped]=ffc_plan_dicom_batch({unknown,svs},tempdir);
    verifyEmpty(tc,p); verifyNumElements(tc,skipped,2);
end
function testMissingStudyDateDoesNotMergeScans(tc)
    a=scanMeta('A','Alpha','');
    p=ffc_plan_dicom_batch({a,a},tempdir);
    verifyNotEqual(tc,p(1).options.StudyInstanceUID,p(2).options.StudyInstanceUID);
end
function testWriterUsesPerScanSnapshot(tc)
    assumeTrue(tc,exist('dicominfo','file')==2);
    % A small struct fixture exercises the writer without reconstruction.
    a=scanMeta('A','Alpha','20261110'); b=scanMeta('B','Beta','20261110');
    root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
    p=ffc_plan_dicom_batch({a,b},root);
    ffc_export_enhanced_mr(a,p(1).filename,p(1).options);
    ffc_export_enhanced_mr(b,p(2).filename,p(2).options);
    first=dicominfo(p(1).filename); second=dicominfo(p(2).filename);
    verifyEqual(tc,first.PatientID,'A'); verifyEqual(tc,second.PatientID,'B');
    verifyNotEqual(tc,first.StudyInstanceUID,second.StudyInstanceUID);
    verifyEqual(tc,first.NumberOfFrames,second.NumberOfFrames);
end
function scan = scanMeta(id,name,date)
    scan=struct('magimage',ones(4,4,2,'single'),'is_svs',false, ...
        'source_metadata',struct('PatientID',id,'PatientName',name,'Session',date, ...
        'StudyID','FCI','StudyDescription','FCI','SerieName','Same series','SerieID','1'), ...
        'param',struct(),'sequence','GRADIENT_ECHO','displayname','Same series', ...
        'filelocation','','orientation',[1 0 0 0 -1 0], ...
        'TwoDimensional',1,'fov',40,'fov_phase',40,'fov_3d',20,'thk',10, ...
        'slices',2,'samples',4,'views',4,'fieldpoints',200,'timepoints',100, ...
        'fft_size',4,'fft_size_3d',2,'TE',10,'FlipAngle',90);
end

function testParallelBatchMatchesSerialPixels(tc)
    assumeTrue(tc,license('test','Distrib_Computing_Toolbox'));
    assumeTrue(tc,exist('dicomread','file')==2);
    pool=gcp('nocreate'); assumeFalse(tc,isa(pool,'parallel.ThreadPool'));
    root=tempname; mkdir(root);
    had=ispref('FFCProcessing','ParallelDicomExport');
    old=[]; if had, old=getpref('FFCProcessing','ParallelDicomExport'); end
    cleanup=onCleanup(@restore); %#ok<NASGU>
    scans=cell(1,2);
    for k=1:2
        scan=makeScan(true); scan=buildimages(ffc_normalize_viewer_scan(scan),false);
        m=scanMeta(sprintf('P%d',k),sprintf('Person%d',k),'20261110');
        scan.source_metadata=m.source_metadata; scan.filelocation=''; scans{k}=scan;
    end
    setpref('FFCProcessing','ParallelDicomExport',false);
    serial=ffc_export_dicom_batch(scans,root);
    setpref('FFCProcessing','ParallelDicomExport',true);
    parallel=ffc_export_dicom_batch(scans,root);
    verifyEmpty(tc,serial.failed); verifyEmpty(tc,parallel.failed);
    verifyEqual(tc,numel(parallel.results),2); verifyGreaterThan(tc,parallel.workers,0);
    for k=1:2
        a=serial.results{find(cellfun(@(x) x.index==k,serial.results),1)};
        b=parallel.results{find(cellfun(@(x) x.index==k,parallel.results),1)};
        verifyEqual(tc,dicomread(a.info.filename),dicomread(b.info.filename));
        verifyEqual(tc,a.info.PatientID,b.info.PatientID);
        verifyEqual(tc,a.info.StudyInstanceUID,b.info.StudyInstanceUID);
    end
    verifyNotEqual(tc,serial.folder,parallel.folder);
    verifyTrue(tc,all(cellfun(@ffc_reconstruction_valid,scans)));
    function restore()
        if had, setpref('FFCProcessing','ParallelDicomExport',old);
        elseif ispref('FFCProcessing','ParallelDicomExport'), rmpref('FFCProcessing','ParallelDicomExport'); end
        if isfolder(root), rmdir(root,'s'); end
    end
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
