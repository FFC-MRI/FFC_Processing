function exportInfo = ffc_export_enhanced_mr(scan, filename, options)
%FFC_EXPORT_ENHANCED_MR Export an FCI reconstruction as Enhanced MR DICOM.
%
%   exportInfo = ffc_export_enhanced_mr(scan, filename, options)
%
% Writes one standards-oriented Enhanced MR Image Storage SOP Instance
% (1.2.840.10008.5.1.4.1.1.4.1) containing the complete reconstructed
% magnitude data set. Frames are ordered by evolution field, evolution time,
% then slice. Evolution field is represented as the DICOM Stack dimension,
% evolution time as TemporalPositionIndex, and slice as InStackPositionNumber.
% The exact FCI evolution field and evolution time are retained in FrameComments
% for every frame.
%
% Geometry is supplied in DICOM LPS by ffc_scan_spatial_geometry. Pixel
% storage reverses both native in-plane dimensions and retains slice order;
% export direction vectors are derived from that exact array transformation.
% No user GUI display transform is baked into the export.
%
% This writer deliberately uses Explicit VR Little Endian and implements the
% small subset of DICOM encoding needed for this Enhanced MR object directly,
% avoiding vendor-private dependencies.

    if nargin < 2 || isempty(filename)
        error('ffc_export_enhanced_mr:FilenameRequired', ...
            'An output DICOM filename is required.');
    end
    if nargin < 3 || isempty(options), options = struct(); end
    if isempty(scan) || ~hasMember(scan,'magimage') || isempty(scan.magimage)
        error('ffc_export_enhanced_mr:NoMagnitudeImage', ...
            'The selected scan does not contain reconstructed magnitude images.');
    end

    if isobject(scan) && isprop(scan, 'reconstruction_key') && ...
            ~ffc_reconstruction_valid(scan)
        error('ffc_export_enhanced_mr:StaleReconstruction', ...
            'Rebuild this scan with TEST25 before export (scan = buildimages(scan)).');
    end

    imageData = double(scan.magimage);
    imageData(~isfinite(imageData)) = 0;
    imageData(imageData < 0) = 0;

    sz = size(imageData);
    sz(end+1:5) = 1;
    nRead = sz(1);
    nPhase = sz(2);
    nSlices = sz(3);
    nTimes = sz(4);
    nFields = sz(5);
    nFrames = nSlices * nTimes * nFields;
    if nFrames < 1
        error('ffc_export_enhanced_mr:NoFrames','No image frames are available.');
    end

    geometry = ffc_scan_spatial_geometry(scan);
    arrayBasis = double(geometry.dimensionBasisPatient);
    spacing = double(geometry.spacingMm(:));
    centre = double(geometry.centerPatientMm(:));
    if size(arrayBasis,1) ~= 3 || size(arrayBasis,2) < 3 || numel(spacing) < 3
        error('ffc_export_enhanced_mr:Geometry', ...
            'The scan does not contain valid patient-space geometry.');
    end

    % Native array is [read, phase, slice]. The written DICOM frame is
    % flip(flip(native,1),2).' (implicit transpose through column-major IO).
    % Increasing DICOM columns/rows therefore DECREASE native dimensions 1/2.
    rows = nPhase;
    columns = nRead;
    rowSpacing = spacing(2);
    columnSpacing = spacing(1);
    sliceSpacing = spacing(3);
    readDirection = unitVector(-arrayBasis(:,1));
    phaseDirection = unitVector(-arrayBasis(:,2));
    sliceDirection = unitVector(arrayBasis(:,3));
    iop = [readDirection(:).' phaseDirection(:).'];

    fieldValues = vectorForCount(memberValue(scan,'fieldpoints',[]),nFields,NaN);
    timeValuesMs = timeMatrixForExport(memberValue(scan,'timepoints',[]),nTimes,nFields);

    maxPixel = max(imageData(:));
    if ~(isfinite(maxPixel) && maxPixel > 0)
        maxPixel = 1;
    end
    storedMax = 65535;
    storedScale = storedMax / maxPixel;
    physicalSlope = 1 / storedScale;

    sopClassUID = '1.2.840.10008.5.1.4.1.1.4.1';
    transferSyntaxUID = '1.2.840.10008.1.2.1';
    if isfield(options,'ResolvedMetadata') && isstruct(options.ResolvedMetadata)
        % Batch planner snapshots each scan independently on the client.
        rawMeta = options.ResolvedMetadata;
    else
        rawMeta = ffc_extract_dicom_metadata(scan);
    end
    if isfield(options,'PatientMetadata') && isstruct(options.PatientMetadata)
        rawMeta = mergeMetadata(rawMeta,options.PatientMetadata);
    end

    sopInstanceUID = optionText(options,'SOPInstanceUID','');
    if isempty(sopInstanceUID), sopInstanceUID = ffc_make_dicom_uid(); end
    studyUID = optionText(options,'StudyInstanceUID','');
    if isempty(studyUID)
        if ~isempty(rawMeta.StudyKey)
            studyUID = ffc_make_dicom_uid(['FCI-STUDY|' rawMeta.StudyKey]);
        else
            studyUID = ffc_make_dicom_uid();
        end
    end
    seriesUID = optionText(options,'SeriesInstanceUID','');
    if isempty(seriesUID), seriesUID = ffc_make_dicom_uid(); end
    frameOfReferenceUID = optionText(options,'FrameOfReferenceUID','');
    if isempty(frameOfReferenceUID)
        if ~isempty(rawMeta.StudyKey)
            frameOfReferenceUID = ffc_make_dicom_uid(['FCI-FOR-LPS25|' rawMeta.StudyKey]);
        else
            frameOfReferenceUID = ffc_make_dicom_uid();
        end
    end
    dimensionOrganizationUID = ffc_make_dicom_uid();
    implementationUID = ffc_make_dicom_uid('ABDN-FCI-ENHANCED-MR-WRITER');

    patientName = safeText(rawMeta.PatientName,64);
    patientID = safeText(rawMeta.PatientID,64);
    patientBirthDate = dicomDate(rawMeta.PatientBirthDate);
    patientSex = upper(rawMeta.PatientSex);
    if ~any(strcmp(patientSex,{'M','F','O'})), patientSex = ''; end

    seriesDescription = safeText(rawMeta.SeriesDescription,64);
    if isempty(rawMeta.SeriesDescription)
        seriesDescription = ['FCI Enhanced MR - ' safeText(memberValue(scan,'displayname','reconstruction'),40)];
    end
    protocolName = safeText(memberValue(scan,'sequence',seriesDescription),64);
    if nSlices > 1
        imageFlavor = 'VOLUME';
    else
        imageFlavor = 'STATIC';
    end
    [nowDate,nowTime,nowDateTime] = currentDicomDateTime();
    studyDate = rawMeta.StudyDate;
    studyTime = rawMeta.StudyTime;
    [acqDate,acqTime,acqDateTime] = normalizeAcquisitionDateTime(rawMeta.AcquisitionDateTime);
    if isempty(studyDate), studyDate = acqDate; end
    if isempty(studyDate), studyDate = nowDate; end
    if isempty(studyTime), studyTime = acqTime; end
    if isempty(studyTime), studyTime = nowTime; end
    if isempty(acqDate), acqDate = studyDate; end
    if isempty(acqTime), acqTime = studyTime; end
    if isempty(acqDateTime), acqDateTime = [acqDate acqTime]; end
    if isempty(acqDateTime), acqDateTime = nowDateTime; end
    seriesDate = acqDate; seriesTime = acqTime;
    contentDate = acqDate; contentTime = acqTime;

    teMs = finiteScalar(memberValue(scan,'TE',NaN),NaN);
    flipAngle = finiteScalar(memberValue(scan,'FlipAngle',NaN),NaN);
    sliceThickness = finiteScalar(memberValue(scan,'thk',NaN),sliceSpacing);
    if ~(isfinite(sliceThickness) && sliceThickness > 0), sliceThickness = sliceSpacing; end

    fid = fopen(filename,'w','ieee-le');
    if fid < 0
        error('ffc_export_enhanced_mr:OpenFailed','Could not open %s for writing.',filename);
    end
    cleanup = onCleanup(@() fcloseIfOpen(fid)); %#ok<NASGU>

    % Part 10 preamble + File Meta Information.
    fwrite(fid,zeros(128,1,'uint8'),'uint8');
    fwrite(fid,uint8('DICM'),'uint8');
    groupLengthValuePos = writeULPlaceholder(fid,hex2dec('0002'),hex2dec('0000'));
    metaStart = ftell(fid);
    writeElement(fid,hex2dec('0002'),hex2dec('0001'),'OB',uint8([0 1]));
    writeElement(fid,hex2dec('0002'),hex2dec('0002'),'UI',sopClassUID);
    writeElement(fid,hex2dec('0002'),hex2dec('0003'),'UI',sopInstanceUID);
    writeElement(fid,hex2dec('0002'),hex2dec('0010'),'UI',transferSyntaxUID);
    writeElement(fid,hex2dec('0002'),hex2dec('0012'),'UI',implementationUID);
    writeElement(fid,hex2dec('0002'),hex2dec('0013'),'SH','ABDNFCI_TEST25');
    metaEnd = ftell(fid);
    currentPos = ftell(fid);
    fseek(fid,groupLengthValuePos,'bof');
    fwrite(fid,uint32(metaEnd-metaStart),'uint32');
    fseek(fid,currentPos,'bof');

    % Patient / Study / Series / Equipment.
    writeElement(fid,hex2dec('0008'),hex2dec('0005'),'CS','ISO_IR 100');
    writeElement(fid,hex2dec('0008'),hex2dec('0008'),'CS',{'DERIVED','PRIMARY',imageFlavor,'NONE'});
    writeElement(fid,hex2dec('0008'),hex2dec('0016'),'UI',sopClassUID);
    writeElement(fid,hex2dec('0008'),hex2dec('0018'),'UI',sopInstanceUID);
    writeElement(fid,hex2dec('0008'),hex2dec('0020'),'DA',studyDate);
    writeElement(fid,hex2dec('0008'),hex2dec('0021'),'DA',seriesDate);
    writeElement(fid,hex2dec('0008'),hex2dec('0022'),'DA',acqDate);
    writeElement(fid,hex2dec('0008'),hex2dec('0023'),'DA',contentDate);
    writeElement(fid,hex2dec('0008'),hex2dec('002A'),'DT',acqDateTime);
    writeElement(fid,hex2dec('0008'),hex2dec('0030'),'TM',studyTime);
    writeElement(fid,hex2dec('0008'),hex2dec('0031'),'TM',seriesTime);
    writeElement(fid,hex2dec('0008'),hex2dec('0032'),'TM',acqTime);
    writeElement(fid,hex2dec('0008'),hex2dec('0033'),'TM',contentTime);
    writeElement(fid,hex2dec('0008'),hex2dec('0050'),'SH','');
    writeElement(fid,hex2dec('0008'),hex2dec('0060'),'CS','MR');
    writeElement(fid,hex2dec('0008'),hex2dec('0070'),'LO','University of Aberdeen');
    writeElement(fid,hex2dec('0008'),hex2dec('0080'),'LO','Aberdeen Biomedical Imaging Centre');
    writeElement(fid,hex2dec('0008'),hex2dec('0090'),'PN','');
    writeElement(fid,hex2dec('0008'),hex2dec('1030'),'LO',safeText(rawMeta.StudyDescription,64));
    writeElement(fid,hex2dec('0008'),hex2dec('103E'),'LO',seriesDescription);
    writeElement(fid,hex2dec('0008'),hex2dec('1090'),'LO','FCI MRI');

    writeElement(fid,hex2dec('0010'),hex2dec('0010'),'PN',patientName);
    writeElement(fid,hex2dec('0010'),hex2dec('0020'),'LO',patientID);
    writeElement(fid,hex2dec('0010'),hex2dec('0030'),'DA',patientBirthDate);
    writeElement(fid,hex2dec('0010'),hex2dec('0040'),'CS',patientSex);
    patientWeight = str2double(rawMeta.PatientWeight);
    if isscalar(patientWeight) && isfinite(patientWeight) && patientWeight > 0
        writeElement(fid,hex2dec('0010'),hex2dec('1030'),'DS',patientWeight);
    end

    writeElement(fid,hex2dec('0018'),hex2dec('0015'),'CS','WHOLEBODY'); % generic fallback anatomy
    writeElement(fid,hex2dec('0018'),hex2dec('0087'),'DS',0.2); % readout field (T)
    writeElement(fid,hex2dec('0018'),hex2dec('1000'),'LO','FCI-RESEARCH');
    writeElement(fid,hex2dec('0018'),hex2dec('1020'),'LO','FFC Processing TEST25');
    writeElement(fid,hex2dec('0018'),hex2dec('1030'),'LO',protocolName);
    writeElement(fid,hex2dec('0018'),hex2dec('5100'),'CS',rawMeta.PatientPosition);

    writeElement(fid,hex2dec('0020'),hex2dec('000D'),'UI',studyUID);
    writeElement(fid,hex2dec('0020'),hex2dec('000E'),'UI',seriesUID);
    studyID = safeText(rawMeta.StudyID,16);
    if isempty(studyID), studyID = safeText(rawMeta.StudyDescription,16); end
    rawSeriesNumber = str2double(rawMeta.SeriesID);
    if ~(isscalar(rawSeriesNumber) && isfinite(rawSeriesNumber) && rawSeriesNumber >= 0), rawSeriesNumber = 9001; end
    writeElement(fid,hex2dec('0020'),hex2dec('0010'),'SH',studyID);
    writeElement(fid,hex2dec('0020'),hex2dec('0011'),'IS',optionNumber(options,'SeriesNumber',rawSeriesNumber));
    writeElement(fid,hex2dec('0020'),hex2dec('0012'),'IS',1);
    writeElement(fid,hex2dec('0020'),hex2dec('0013'),'IS',1);
    writeElement(fid,hex2dec('0020'),hex2dec('0105'),'IS',nTimes);
    writeElement(fid,hex2dec('0020'),hex2dec('0052'),'UI',frameOfReferenceUID);
    writeElement(fid,hex2dec('0020'),hex2dec('1040'),'LO','');

    % Enhanced MR image pixel module.
    writeElement(fid,hex2dec('0028'),hex2dec('0002'),'US',1);
    writeElement(fid,hex2dec('0028'),hex2dec('0004'),'CS','MONOCHROME2');
    writeElement(fid,hex2dec('0028'),hex2dec('0008'),'IS',nFrames);
    writeElement(fid,hex2dec('0028'),hex2dec('0010'),'US',rows);
    writeElement(fid,hex2dec('0028'),hex2dec('0011'),'US',columns);
    writeElement(fid,hex2dec('0028'),hex2dec('0100'),'US',16);
    writeElement(fid,hex2dec('0028'),hex2dec('0101'),'US',16);
    writeElement(fid,hex2dec('0028'),hex2dec('0102'),'US',15);
    writeElement(fid,hex2dec('0028'),hex2dec('0103'),'US',0);
    writeElement(fid,hex2dec('0028'),hex2dec('0301'),'CS','NO');
    writeElement(fid,hex2dec('0028'),hex2dec('2110'),'CS','00');
    writeElement(fid,hex2dec('2050'),hex2dec('0020'),'CS','IDENTITY');

    % Common Enhanced MR image description attributes.
    writeElement(fid,hex2dec('0018'),hex2dec('9004'),'CS','RESEARCH'); % Content Qualification
    writeElement(fid,hex2dec('0018'),hex2dec('9100'),'CS','1H'); % Resonant Nucleus
    writeElement(fid,hex2dec('0008'),hex2dec('9205'),'CS','MONOCHROME'); % Pixel Presentation
    writeElement(fid,hex2dec('0008'),hex2dec('9206'),'CS','VOLUME');     % Volumetric Properties
    writeElement(fid,hex2dec('0008'),hex2dec('9207'),'CS','NONE');       % Volume Based Calculation Technique
    writeElement(fid,hex2dec('0008'),hex2dec('9208'),'CS','MAGNITUDE');  % Complex Image Component
    writeElement(fid,hex2dec('0008'),hex2dec('9209'),'CS','UNKNOWN');    % Acquisition Contrast

    % Multi-frame dimensions: evolution field (stack), evolution time and slice.
    % The Dimension Index Values themselves are ordinal, as required by DICOM;
    % exact physical field/time values are retained in Frame Comments.
    beginSequence(fid,hex2dec('0020'),hex2dec('9221')); % Dimension Organization Sequence
        beginItem(fid);
            writeElement(fid,hex2dec('0020'),hex2dec('9164'),'UI',dimensionOrganizationUID);
        endItem(fid);
    endSequence(fid);

    beginSequence(fid,hex2dec('0020'),hex2dec('9222')); % Dimension Index Sequence
        beginItem(fid);
            writeElement(fid,hex2dec('0020'),hex2dec('9164'),'UI',dimensionOrganizationUID);
            writeElement(fid,hex2dec('0020'),hex2dec('9165'),'AT',[hex2dec('0020') hex2dec('9056')]);
            writeElement(fid,hex2dec('0020'),hex2dec('9167'),'AT',[hex2dec('0020') hex2dec('9111')]);
            writeElement(fid,hex2dec('0020'),hex2dec('9421'),'LO','FCI evolution field');
        endItem(fid);
        beginItem(fid);
            writeElement(fid,hex2dec('0020'),hex2dec('9164'),'UI',dimensionOrganizationUID);
            writeElement(fid,hex2dec('0020'),hex2dec('9165'),'AT',[hex2dec('0020') hex2dec('9128')]);
            writeElement(fid,hex2dec('0020'),hex2dec('9167'),'AT',[hex2dec('0020') hex2dec('9111')]);
            writeElement(fid,hex2dec('0020'),hex2dec('9421'),'LO','FCI evolution time');
        endItem(fid);
        beginItem(fid);
            writeElement(fid,hex2dec('0020'),hex2dec('9164'),'UI',dimensionOrganizationUID);
            writeElement(fid,hex2dec('0020'),hex2dec('9165'),'AT',[hex2dec('0020') hex2dec('9057')]);
            writeElement(fid,hex2dec('0020'),hex2dec('9167'),'AT',[hex2dec('0020') hex2dec('9111')]);
            writeElement(fid,hex2dec('0020'),hex2dec('9421'),'LO','Slice');
        endItem(fid);
    endSequence(fid);

    % Acquisition Context module is mandatory; an empty sequence is valid
    % when there is no coded acquisition-context content to add.
    writeEmptySequence(fid,hex2dec('0040'),hex2dec('0555'));

    % Shared Functional Groups Sequence.
    beginSequence(fid,hex2dec('5200'),hex2dec('9229'));
        beginItem(fid);
            beginSequence(fid,hex2dec('0028'),hex2dec('9110')); % Pixel Measures
                beginItem(fid);
                    writeElement(fid,hex2dec('0028'),hex2dec('0030'),'DS',[rowSpacing columnSpacing]);
                    writeElement(fid,hex2dec('0018'),hex2dec('0050'),'DS',sliceThickness);
                    writeElement(fid,hex2dec('0018'),hex2dec('0088'),'DS',sliceSpacing);
                endItem(fid);
            endSequence(fid);

            beginSequence(fid,hex2dec('0020'),hex2dec('9116')); % Plane Orientation (Patient)
                beginItem(fid);
                    writeElement(fid,hex2dec('0020'),hex2dec('0037'),'DS',iop);
                endItem(fid);
            endSequence(fid);

            % Frame Anatomy is mandatory for Enhanced MR. When no more
            % specific anatomy is available from the FCI acquisition metadata, use
            % the DICOM/SNOMED 'Entire body' concept as a standards-valid generic
            % fallback. Frame Laterality U denotes an unpaired region.
            beginSequence(fid,hex2dec('0020'),hex2dec('9071')); % Frame Anatomy
                beginItem(fid);
                    writeElement(fid,hex2dec('0020'),hex2dec('9072'),'CS','U');
                    beginSequence(fid,hex2dec('0008'),hex2dec('2218')); % Anatomic Region
                        beginItem(fid);
                            writeElement(fid,hex2dec('0008'),hex2dec('0100'),'SH','38266002');
                            writeElement(fid,hex2dec('0008'),hex2dec('0102'),'SH','SCT');
                            writeElement(fid,hex2dec('0008'),hex2dec('0104'),'LO','Entire body');
                        endItem(fid);
                    endSequence(fid);
                endItem(fid);
            endSequence(fid);

            beginSequence(fid,hex2dec('0018'),hex2dec('9226')); % MR Image Frame Type
                beginItem(fid);
                    writeElement(fid,hex2dec('0008'),hex2dec('9007'),'CS',{'DERIVED','PRIMARY',imageFlavor,'NONE'});
                    writeElement(fid,hex2dec('0008'),hex2dec('9205'),'CS','MONOCHROME');
                    writeElement(fid,hex2dec('0008'),hex2dec('9206'),'CS','VOLUME');
                    writeElement(fid,hex2dec('0008'),hex2dec('9207'),'CS','NONE');
                    writeElement(fid,hex2dec('0008'),hex2dec('9208'),'CS','MAGNITUDE');
                    writeElement(fid,hex2dec('0008'),hex2dec('9209'),'CS','UNKNOWN');
                endItem(fid);
            endSequence(fid);

            if isfinite(teMs) && teMs >= 0
                beginSequence(fid,hex2dec('0018'),hex2dec('9114')); % MR Echo
                    beginItem(fid);
                        writeElement(fid,hex2dec('0018'),hex2dec('9082'),'FD',teMs); % Effective Echo Time
                    endItem(fid);
                endSequence(fid);
            end

            if isfinite(flipAngle) && flipAngle >= 0
                beginSequence(fid,hex2dec('0018'),hex2dec('9112')); % MR Timing and Related Parameters
                    beginItem(fid);
                        writeElement(fid,hex2dec('0018'),hex2dec('1314'),'DS',flipAngle);
                    endItem(fid);
                endSequence(fid);
            end

            beginSequence(fid,hex2dec('0028'),hex2dec('9145')); % Pixel Value Transformation
                beginItem(fid);
                    writeElement(fid,hex2dec('0028'),hex2dec('1052'),'DS',0);
                    writeElement(fid,hex2dec('0028'),hex2dec('1053'),'DS',physicalSlope);
                    writeElement(fid,hex2dec('0028'),hex2dec('1054'),'LO','US');
                endItem(fid);
            endSequence(fid);

            beginSequence(fid,hex2dec('0028'),hex2dec('9132')); % Frame VOI LUT
                beginItem(fid);
                    writeElement(fid,hex2dec('0028'),hex2dec('1050'),'DS',maxPixel*0.35);
                    writeElement(fid,hex2dec('0028'),hex2dec('1051'),'DS',maxPixel*0.70);
                    writeElement(fid,hex2dec('0028'),hex2dec('1055'),'LO','FCI magnitude');
                endItem(fid);
            endSequence(fid);
        endItem(fid);
    endSequence(fid);

    % Per-frame geometry and FCI condition metadata.
    beginSequence(fid,hex2dec('5200'),hex2dec('9230'));
    frameNumber = 0;
    for fieldIndex = 1:nFields
        for timeIndex = 1:nTimes
            temporalIndex = timeIndex;
            stackID = sprintf('F%03d',fieldIndex);
            for sliceIndex = 1:nSlices
                frameNumber = frameNumber + 1;
                sliceOffset = (sliceIndex - (nSlices+1)/2) * sliceSpacing;
                sliceCenter = centre + sliceOffset*sliceDirection;
                imagePosition = sliceCenter ...
                    - ((columns-1)/2)*columnSpacing*readDirection ...
                    - ((rows-1)/2)*rowSpacing*phaseDirection;

                field_mT = fieldValues(min(fieldIndex,numel(fieldValues)));
                time_ms = timeValuesMs(timeIndex,min(fieldIndex,size(timeValuesMs,2)));
                comments = sprintf('FCI evolution field=%s mT; evolution time=%s ms', ...
                    numberText(field_mT),numberText(time_ms));

                beginItem(fid);
                    beginSequence(fid,hex2dec('0020'),hex2dec('9111')); % Frame Content
                        beginItem(fid);
                            writeElement(fid,hex2dec('0020'),hex2dec('9157'),'UL',[fieldIndex temporalIndex sliceIndex]);
                            writeElement(fid,hex2dec('0020'),hex2dec('9056'),'SH',stackID);
                            writeElement(fid,hex2dec('0020'),hex2dec('9057'),'UL',sliceIndex);
                            writeElement(fid,hex2dec('0020'),hex2dec('9128'),'UL',temporalIndex);
                            writeElement(fid,hex2dec('0020'),hex2dec('9156'),'US',frameNumber);
                            writeElement(fid,hex2dec('0020'),hex2dec('9158'),'LT',comments);
                        endItem(fid);
                    endSequence(fid);
                    beginSequence(fid,hex2dec('0020'),hex2dec('9113')); % Plane Position (Patient)
                        beginItem(fid);
                            writeElement(fid,hex2dec('0020'),hex2dec('0032'),'DS',imagePosition(:).');
                        endItem(fid);
                    endSequence(fid);
                endItem(fid);
            end
        end
    end
    endSequence(fid);

    % Pixel Data (OW). Because the DICOM image for each frame is the
    % transpose of the native MATLAB matrix, writing frame(:) in MATLAB's
    % column-major order is exactly DICOM row-major order for that transpose.
    pixelBytes = double(nRead) * double(nPhase) * double(nFrames) * 2;
    if pixelBytes >= double(2^32-2)
        error('ffc_export_enhanced_mr:TooLarge', ...
            'Pixel Data exceeds the Explicit-VR single-element size supported by this exporter.');
    end
    writeElementHeader(fid,hex2dec('7FE0'),hex2dec('0010'),'OW',uint32(pixelBytes));
    for fieldIndex = 1:nFields
        for timeIndex = 1:nTimes
            for sliceIndex = 1:nSlices
                % Geometry above accounts for these two pixel reversals.
                sourceSliceIndex = sliceIndex;
                frame = imageData(:,:,sourceSliceIndex,timeIndex,fieldIndex);
                frame = flip(flip(frame,1),2);
                stored = uint16(min(storedMax,max(0,round(frame*storedScale))));
                fwrite(fid,stored(:),'uint16');
            end
        end
    end

    exportInfo = struct();
    exportInfo.filename = filename;
    exportInfo.SOPClassUID = sopClassUID;
    exportInfo.SOPInstanceUID = sopInstanceUID;
    exportInfo.StudyInstanceUID = studyUID;
    exportInfo.SeriesInstanceUID = seriesUID;
    exportInfo.FrameOfReferenceUID = frameOfReferenceUID;
    exportInfo.NumberOfFrames = nFrames;
    exportInfo.NativeMatrix = [nRead nPhase nSlices nTimes nFields];
    exportInfo.DICOMRowsColumns = [rows columns];
    exportInfo.PixelSpacing = [rowSpacing columnSpacing];
    exportInfo.ImageOrientationPatient = iop;
    exportInfo.CoordinateConvention = 'TEST25: shared voxel geometry in DICOM LPS; Cameleon Y/Z conversion; both in-plane indices reversed; native slice order';
    exportInfo.GlobalMagnitudeToStoredScale = storedScale;
    exportInfo.StoredToMagnitudeSlope = physicalSlope;
    exportInfo.FrameOrder = 'field -> time -> slice';
    exportInfo.FCIConditionEncoding = ...
        'Stack ID=evolution field index; TemporalPositionIndex=evolution time index; exact values in FrameComments';
    exportInfo.PatientName = patientName;
    exportInfo.PatientID = patientID;
    exportInfo.PatientBirthDate = patientBirthDate;
    exportInfo.PatientSex = patientSex;
    exportInfo.StudyID = studyID;
    exportInfo.StudyDate = studyDate;
    exportInfo.AcquisitionDateTime = acqDateTime;
    exportInfo.StudyDescription = rawMeta.StudyDescription;
    exportInfo.SeriesDescription = seriesDescription;
    exportInfo.SeriesNumber = optionNumber(options,'SeriesNumber',9001);
    exportInfo.MetadataSource = rawMeta.MetadataSource;
end

function merged = mergeMetadata(base,override)
    merged = base;
    fields = fieldnames(override);
    for k = 1:numel(fields)
        v = override.(fields{k});
        if ischar(v) || isstring(v)
            if strlength(string(v)) == 0, continue; end
        elseif isempty(v)
            continue
        end
        merged.(fields{k}) = v;
    end
end

function value = optionText(options,name,defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options,name)
        c = options.(name);
        if isstring(c), c = char(c(1)); end
        if ischar(c) && ~isempty(strtrim(c)), value = strtrim(c); end
    end
end

function value = optionNumber(options,name,defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options,name)
        c = double(options.(name));
        if isscalar(c) && isfinite(c), value = c; end
    end
end

function [d,t,dt] = normalizeAcquisitionDateTime(raw)
    d = ''; t = ''; dt = '';
    raw = char(string(raw));
    m = regexp(raw,'(\d{4})[-/]?(\d{2})[-/]?(\d{2})','tokens','once');
    if ~isempty(m), d = [m{1} m{2} m{3}]; end
    m = regexp(raw,'[T ](\d{2}):?(\d{2}):?(\d{2})(?:\.(\d+))?','tokens','once');
    if ~isempty(m)
        t = [m{1} m{2} m{3}];
        if numel(m) >= 4 && ~isempty(m{4}), t = [t '.' m{4}]; end
    end
    if ~isempty(d) && ~isempty(t), dt = [d t]; end
end

function value = memberValue(obj,name,defaultValue)
    value = defaultValue;
    try
        if (isobject(obj) && isprop(obj,name)) || (isstruct(obj) && isfield(obj,name))
            candidate = obj.(name);
            if ~isempty(candidate), value = candidate; end
        end
    catch
    end
end

function tf = hasMember(obj,name)
    tf = false;
    try
        tf = (isobject(obj) && isprop(obj,name)) || (isstruct(obj) && isfield(obj,name));
    catch
    end
end

function value = finiteScalar(candidate,defaultValue)
    value = defaultValue;
    if isnumeric(candidate) || islogical(candidate)
        c = double(candidate(:));
        c = c(isfinite(c));
        if ~isempty(c), value = c(1); end
    end
end

function v = unitVector(v)
    v = double(v(:));
    n = norm(v);
    if ~(isfinite(n) && n > 0), error('Invalid geometry direction vector.'); end
    v = v./n;
end

function values = vectorForCount(candidate,count,defaultValue)
    if isnumeric(candidate) || islogical(candidate)
        c = double(candidate(:).');
    else
        c = [];
    end
    if isempty(c)
        values = repmat(defaultValue,1,count);
    elseif numel(c) == 1
        values = repmat(c,1,count);
    elseif numel(c) >= count
        values = c(1:count);
    else
        values = [c repmat(c(end),1,count-numel(c))];
    end
end

function t = timeMatrixForExport(candidate,nTimes,nFields)
    if ~(isnumeric(candidate) || islogical(candidate)) || isempty(candidate)
        t = nan(nTimes,nFields);
        return
    end
    c = double(candidate);
    if isvector(c)
        c = c(:);
        if numel(c) >= nTimes
            t = repmat(c(1:nTimes),1,nFields);
        else
            c(end+1:nTimes,1) = c(end);
            t = repmat(c,1,nFields);
        end
        return
    end
    if size(c,1) == nFields && size(c,2) == nTimes
        c = c.';
    end
    t = nan(nTimes,nFields);
    nr = min(nTimes,size(c,1));
    nc = min(nFields,size(c,2));
    t(1:nr,1:nc) = c(1:nr,1:nc);
    for f = 1:nFields
        finiteCol = t(isfinite(t(:,f)),f);
        if isempty(finiteCol)
            if f > 1, t(:,f) = t(:,f-1); end
        elseif numel(finiteCol) < nTimes
            last = finiteCol(end);
            t(~isfinite(t(:,f)),f) = last;
        end
    end
end

function value = numericParameter(scan,name,defaultValue)
    value = defaultValue;
    params = memberValue(scan,'param',struct());
    if isstruct(params) && isfield(params,name)
        candidate = params.(name);
        if isnumeric(candidate) || islogical(candidate)
            candidate = double(candidate(:));
            candidate = candidate(isfinite(candidate));
            if ~isempty(candidate)
                value = candidate(1);
                return
            end
        end
        try
            candidate = str2double(string(candidate));
            candidate = candidate(isfinite(candidate));
            if ~isempty(candidate)
                value = double(candidate(1));
            end
        catch
        end
    end
end

function text = textParameter(scan,names,defaultValue)
    text = defaultValue;
    candidates = {};
    for k = 1:numel(names)
        candidates{end+1} = memberValue(scan,names{k},[]); %#ok<AGROW>
    end
    params = memberValue(scan,'param',struct());
    if isstruct(params)
        for k = 1:numel(names)
            if isfield(params,names{k})
                candidates{end+1} = params.(names{k}); %#ok<AGROW>
            end
        end
    end
    for k = 1:numel(candidates)
        c = candidates{k};
        if iscell(c) && ~isempty(c), c = c{1}; end
        if isstring(c), c = char(c(1)); end
        if ischar(c) && ~isempty(strtrim(c))
            text = strtrim(c);
            return
        end
    end
end

function out = safeText(value,maxLen)
    if iscell(value) && ~isempty(value), value = value{1}; end
    if isstring(value), value = char(value(1)); end
    if ~ischar(value), value = 'FCI'; end
    out = regexprep(value,'[\\\r\n\t]',' ');
    out = strtrim(out);
    if isempty(out), out = 'FCI'; end
    if nargin > 1 && numel(out) > maxLen, out = out(1:maxLen); end
end

function out = dicomDate(value)
    value = safeText(value,32);
    if strcmp(value,'FCI'), value = ''; end
    digitsOnly = regexprep(value,'[^0-9]','');
    if numel(digitsOnly) >= 8, out = digitsOnly(1:8); else, out = ''; end
end

function [d,t,dt] = currentDicomDateTime()
    c = clock;
    d = sprintf('%04d%02d%02d',c(1),c(2),c(3));
    sec = c(6);
    t = sprintf('%02d%02d%06.3f',c(4),c(5),sec);
    dt = [d t];
end

function uid = makeUID()
    try
        if exist('dicomuid','file') == 2
            uid = dicomuid;
            return
        end
    catch
    end
    try
        jUUID = javaMethod('randomUUID','java.util.UUID');
        uuid = char(jUUID.toString());
        hex = strrep(uuid,'-','');
        big = javaObject('java.math.BigInteger',hex,16);
        uid = ['2.25.' char(big.toString())];
    catch
        % Last-resort numeric UID using the 2.25 UUID namespace form. This is
        % not used when Image Processing Toolbox or Java is available.
        c = clock;
        seed = abs(round(sum(c.*[1e10 1e8 1e6 1e4 1e2 1])));
        uid = ['2.25.' sprintf('%.0f',seed) sprintf('%09d',randi(1e9-1))];
    end
end

function s = numberText(v)
    if isfinite(v), s = sprintf('%.10g',v); else, s = 'unknown'; end
end

function fcloseIfOpen(fid)
    try
        if fid > 0, fclose(fid); end
    catch
    end
end

% ----------------------------- DICOM writer -----------------------------
function valuePos = writeULPlaceholder(fid,group,element)
    writeTag(fid,group,element);
    fwrite(fid,uint8('UL'),'uint8');
    fwrite(fid,uint16(4),'uint16');
    valuePos = ftell(fid);
    fwrite(fid,uint32(0),'uint32');
end

function writeElement(fid,group,element,vr,value)
    [bytes,isRawNumeric,valueLength] = encodeValue(vr,value);
    writeElementHeader(fid,group,element,vr,uint32(valueLength));
    if isRawNumeric
        fwrite(fid,bytes.data,bytes.precision);
    else
        fwrite(fid,bytes,'uint8');
    end
end

function writeElementHeader(fid,group,element,vr,valueLength)
    writeTag(fid,group,element);
    fwrite(fid,uint8(vr),'uint8');
    longVR = any(strcmp(vr,{'OB','OD','OF','OL','OV','OW','SQ','UC','UR','UT','UN'}));
    if longVR
        fwrite(fid,uint16(0),'uint16');
        fwrite(fid,uint32(valueLength),'uint32');
    else
        if double(valueLength) > double(intmax('uint16'))
            error('DICOM element %04X,%04X is too long for VR %s.',group,element,vr);
        end
        fwrite(fid,uint16(valueLength),'uint16');
    end
end

function [out,isRawNumeric,valueLength] = encodeValue(vr,value)
    isRawNumeric = false;
    switch vr
        case {'US'}
            data = uint16(value(:));
            out = struct('data',data,'precision','uint16'); isRawNumeric = true; valueLength = 2*numel(data);
        case {'SS'}
            data = int16(value(:));
            out = struct('data',data,'precision','int16'); isRawNumeric = true; valueLength = 2*numel(data);
        case {'UL'}
            data = uint32(value(:));
            out = struct('data',data,'precision','uint32'); isRawNumeric = true; valueLength = 4*numel(data);
        case {'SL'}
            data = int32(value(:));
            out = struct('data',data,'precision','int32'); isRawNumeric = true; valueLength = 4*numel(data);
        case {'FL'}
            data = single(value(:));
            out = struct('data',data,'precision','single'); isRawNumeric = true; valueLength = 4*numel(data);
        case {'FD'}
            data = double(value(:));
            out = struct('data',data,'precision','double'); isRawNumeric = true; valueLength = 8*numel(data);
        case 'AT'
            v = uint16(value);
            if size(v,2) ~= 2, v = reshape(v,2,[]).'; end
            data = reshape(v.',[],1);
            out = struct('data',data,'precision','uint16'); isRawNumeric = true; valueLength = 2*numel(data);
        case 'OB'
            bytes = uint8(value(:));
            if rem(numel(bytes),2), bytes(end+1)=0; end
            out = bytes; valueLength = numel(bytes);
        otherwise
            if isnumeric(value) || islogical(value)
                if strcmp(vr,'DS')
                    parts = arrayfun(@(x) sprintf('%.12g',double(x)),value(:).','UniformOutput',false);
                    text = strjoin(parts,'\');
                elseif strcmp(vr,'IS')
                    parts = arrayfun(@(x) sprintf('%.0f',double(x)),value(:).','UniformOutput',false);
                    text = strjoin(parts,'\');
                else
                    text = num2str(value);
                end
            elseif iscell(value)
                parts = cell(size(value));
                for k = 1:numel(value), parts{k} = char(value{k}); end
                text = strjoin(parts,'\');
            elseif isstring(value)
                text = char(join(value,'\'));
            else
                text = char(value);
            end
            bytes = uint8(text(:).');
            if rem(numel(bytes),2)
                if strcmp(vr,'UI'), bytes(end+1)=uint8(0); else, bytes(end+1)=uint8(' '); end
            end
            out = bytes; valueLength = numel(bytes);
    end
end

function writeTag(fid,group,element)
    fwrite(fid,uint16(group),'uint16');
    fwrite(fid,uint16(element),'uint16');
end

function beginSequence(fid,group,element)
    writeTag(fid,group,element);
    fwrite(fid,uint8('SQ'),'uint8');
    fwrite(fid,uint16(0),'uint16');
    fwrite(fid,uint32(hex2dec('FFFFFFFF')),'uint32');
end

function endSequence(fid)
    writeTag(fid,hex2dec('FFFE'),hex2dec('E0DD'));
    fwrite(fid,uint32(0),'uint32');
end

function writeEmptySequence(fid,group,element)
    writeTag(fid,group,element);
    fwrite(fid,uint8('SQ'),'uint8');
    fwrite(fid,uint16(0),'uint16');
    fwrite(fid,uint32(0),'uint32');
end

function beginItem(fid)
    writeTag(fid,hex2dec('FFFE'),hex2dec('E000'));
    fwrite(fid,uint32(hex2dec('FFFFFFFF')),'uint32');
end

function endItem(fid)
    writeTag(fid,hex2dec('FFFE'),hex2dec('E00D'));
    fwrite(fid,uint32(0),'uint32');
end
