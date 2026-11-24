function tests = test_geometry_test25
% TEST25 geometry and acquisition-offset regressions. Run with runtests.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    root = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.originalPath = path;
    addpath(root); ffc_setup();
end

function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end

function testAxialLPSAndSliceSpacing(testCase)
    scan = fixture();
    g = ffc_scan_spatial_geometry(scan);
    verifyEqual(testCase, g.prescriptionBasisPatient, eye(3), 'AbsTol', 1e-12);
    verifyEqual(testCase, g.dimensionBasisPatient, diag([-1 -1 1]), 'AbsTol', 1e-12);
    verifyEqual(testCase, g.spacingMm, [3.90625 3.90625 12], 'AbsTol', 1e-12);
    firstCentre = g.centerPatientMm - 7*g.spacingMm(3)*g.dimensionBasisPatient(:,3);
    verifyEqual(testCase, firstCentre, [0;0;-84], 'AbsTol', 1e-12);
end

function testAcquisitionHandledOffsetIsNotAppliedTwice(testCase)
    s = rampFixture();
    s.param.OFF_CENTER_FIELD_OF_VIEW_EFF = [-0.025779523882 0 0];
    out = correct_orientation(s);
    ky = (-3:2)/0.060;
    expected = repmat(exp(-2i*pi*ky*0.010), 8, 1);
    verifyEqual(testCase, out.complexkspace, expected, 'AbsTol', 1e-12);
    verifyEqual(testCase, out.fov_offset_info.offsets_m, [0 0.010 0], 'AbsTol', 1e-12);
    again = correct_orientation(out);
    verifyEqual(testCase, again.complexkspace, out.complexkspace);
end

function testMissingEffectiveOffsetsPreservesLegacyRamp(testCase)
    s = rampFixture(); out = correct_orientation(s);
    kx = (-4:3)'/0.080; ky = (-3:2)/0.060;
    expected = exp(-2i*pi*(kx*(-0.025779523882) + ky*0.010));
    verifyEqual(testCase, out.complexkspace, expected, 'AbsTol', 1e-12);
end

function testMalformedEffectiveOffsetsAreRejected(testCase)
    s = rampFixture(); s.param.OFF_CENTER_FIELD_OF_VIEW_EFF = [NaN 0 0];
    verifyError(testCase, @() correct_orientation(s), 'correct_orientation:InvalidEffectiveOffsets');
end

function testEffectivePartitionOffsetOnlyAffects3D(testCase)
    s = rampFixture(); s.complexkspace = complex(ones(8,6,4));
    s.param.OFF_CENTER_FIELD_OF_VIEW_1D = 0;
    s.param.OFF_CENTER_FIELD_OF_VIEW_2D = 0;
    s.param.OFF_CENTER_FIELD_OF_VIEW_3D = 0.012;
    s.param.OFF_CENTER_FIELD_OF_VIEW_EFF = [0 0 0.002];
    twoD = correct_orientation(s);
    verifyEqual(testCase, twoD.complexkspace, s.complexkspace);
    s.TwoDimensional = 0; threeD = correct_orientation(s);
    expected = repmat(reshape(exp(-2i*pi*(-2:1)/0.04*0.010),1,1,4),8,6,1);
    verifyEqual(testCase, threeD.complexkspace, expected, 'AbsTol', 1e-12);
end

function testObliqueCentreAndNativeVoxelMapping(testCase)
    s = fixture(); s.slices = 8; s.magimage = zeros(128,128,8);
    s.orientation = [-0.586165 -0.047265 -0.808812 0.387668 0.860236 -0.331222];
    s.param.IMAGE_ORIENTATION_SUBJECT = s.orientation;
    s.param.SLICE_THICKNESS = 0.007;
    s.param.OFF_CENTER_FIELD_OF_VIEW_X = -0.03845;
    s.param.OFF_CENTER_FIELD_OF_VIEW_Y = 0.08657;
    s.param.OFF_CENTER_FIELD_OF_VIEW_Z = 0.05468;
    g = ffc_scan_spatial_geometry(s);
    verifyEqual(testCase, g.centerPatientMm, [-38.45;-86.57;-54.68], 'AbsTol', 1e-10);
    verifyEqual(testCase, g.spacingMm(3), 9, 'AbsTol', 1e-12);
    r = -g.dimensionBasisPatient(:,1); c = -g.dimensionBasisPatient(:,2);
    n = g.dimensionBasisPatient(:,3);
    ipp = g.centerPatientMm - 3.5*9*n - 63.5*3.90625*(r+c);
    verifyEqual(testCase, ipp, [-11.623275;99.092281;-352.768239], 'AbsTol', 1e-5);
    % Native (read=18,phase=32,slice=3) becomes DICOM (column=111,row=97).
    nativePoint = g.centerPatientMm + g.dimensionBasisPatient * ...
        (g.spacingMm(:).*[18-64.5;32-64.5;3-4.5]);
    frameIPP = ipp + 2*9*n;
    dicomPoint = frameIPP + 110*3.90625*r + 96*3.90625*c;
    verifyEqual(testCase, nativePoint, dicomPoint, 'AbsTol', 1e-10);
end

function testWrittenPixelsAndPositionsAgree(testCase)
    % End-to-end check of the real writer; requires Image Processing Toolbox.
    assumeTrue(testCase, exist('dicomread','file') ~= 0 && exist('dicominfo','file') ~= 0);
    s = fixture(); s.magimage = reshape(1:24, [4 3 2]); s.slices = 2;
    s.fov = 40; s.fov_phase = 30;
    filename = [tempname '.dcm']; cleanup = onCleanup(@() delete(filename)); %#ok<NASGU>
    exported = ffc_export_enhanced_mr(s, filename);
    info = dicominfo(filename); pixels = dicomread(info);
    expected = uint16(round(flip(flip(s.magimage(:,:,1),1),2).' * (65535/24)));
    verifyEqual(testCase, pixels(:,:,1,1), expected);
    verifyEqual(testCase, exported.ImageOrientationPatient, [1 0 0 0 1 0], 'AbsTol', 1e-12);
    position = info.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;
    verifyEqual(testCase, double(position(:)), [-15;-10;-6], 'AbsTol', 1e-10);
end

function s = fixture()
    s = struct('magimage',zeros(128,128,15),'orientation',[1 0 0 0 -1 0], ...
        'fov',500,'fov_phase',500,'fov_3d',178,'thk',10,'slices',15, ...
        'TwoDimensional',1,'is_svs',false);
    s.param = struct('IMAGE_ORIENTATION_SUBJECT',s.orientation, ...
        'SLICE_THICKNESS',0.010,'SPACING_BETWEEN_SLICE',0.002);
end

function s = rampFixture()
    s = struct('complexkspace',complex(ones(8,6)), 'TwoDimensional',1, ...
        'fov_offsets_applied',false,'fov_offset_info',struct(), ...
        'fov',80,'fov_phase',60,'fov_3d',40);
    s.param = struct('OFF_CENTER_FIELD_OF_VIEW_1D',-0.025779523882, ...
        'OFF_CENTER_FIELD_OF_VIEW_2D',0.010,'OFF_CENTER_FIELD_OF_VIEW_3D',0, ...
        'FIELD_OF_VIEW',0.080,'FIELD_OF_VIEW_PHASE',0.060,'FIELD_OF_VIEW_3D',0.040);
end
