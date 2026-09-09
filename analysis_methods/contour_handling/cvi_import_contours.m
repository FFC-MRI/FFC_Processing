function cvi = cvi_import_contours(workspaceFile,dicomRoot)
%CVI_IMPORT_CONTOURS Import cvi42 contours and attach clinical DICOM geometry.
%
%   CVI = CVI_IMPORT_CONTOURS(WORKSPACEFILE,DICOMROOT) reads ventricular
%   contours and converts every matched point into the same subject-coordinate
%   basis used by RS2D IMAGE_ORIENTATION_SUBJECT, in millimetres.  For the
%   supplied cardiac data the RS2D orientation vectors are numerically in the
%   same LPS-like basis as the clinical DICOM direction cosines, so the DICOM
%   X/Y coordinates must NOT be sign-flipped.  TEST6 incorrectly converted
%   them to RAH, producing the large artificial registration rotations.

    cvi = cvi_read_workspace(workspaceFile);
    if ~isfolder(dicomRoot)
        error('cvi_import_contours:DicomFolderNotFound', ...
            'Clinical DICOM folder not found: %s',dicomRoot);
    end
    cvi.dicomRoot = dicomRoot;

    files = dir(fullfile(dicomRoot,'**','*'));
    files = files(~[files.isdir]);
    if isempty(files)
        error('cvi_import_contours:NoDicomFiles','The selected folder is empty.');
    end
    paths = arrayfun(@(x) fullfile(x.folder,x.name),files,'UniformOutput',false);
    uniqueUIDs = unique({cvi.contours.ImageUID},'stable');

    % Find which DICOM image object contains each cvi42 per-frame image UID.
    % The UIDs live in metadata before PixelData in Philips enhanced DICOM,
    % so limiting the scan prevents large cine pixel payloads being read.
    uidMatches = cell(size(uniqueUIDs));
    scanBytes = 2 * 1024 * 1024;
    for f = 1:numel(paths)
        fid = fopen(paths{f},'rb');
        if fid < 0, continue; end
        raw = fread(fid,scanBytes,'*uint8').';
        fclose(fid);
        if isempty(raw), continue; end
        for u = 1:numel(uniqueUIDs)
            if isempty(uidMatches{u}) && ~isempty(strfind(raw,uint8(uniqueUIDs{u}))) %#ok<STREMP>
                uidMatches{u}{end+1} = paths{f}; %#ok<AGROW>
            elseif ~isempty(strfind(raw,uint8(uniqueUIDs{u}))) %#ok<STREMP>
                uidMatches{u}{end+1} = paths{f}; %#ok<AGROW>
            end
        end
    end

    infoCache = containers.Map('KeyType','char','ValueType','any');
    for u = 1:numel(uniqueUIDs)
        candidates = uidMatches{u};
        for f = 1:numel(candidates)
            path = candidates{f};
            if ~isKey(infoCache,path)
                try
                    infoCache(path) = dicominfo(path);
                catch
                    infoCache(path) = [];
                end
            end
        end
    end

    unmatched = false(numel(cvi.contours),1);
    for c = 1:numel(cvi.contours)
        uidIndex = find(strcmp(uniqueUIDs,cvi.contours(c).ImageUID),1);
        candidates = uidMatches{uidIndex};
        selected = '';
        selectedInfo = [];
        for f = 1:numel(candidates)
            info = infoCache(candidates{f});
            if isempty(info) || ~isfield(info,'Rows') || ~isfield(info,'Columns')
                continue;
            end
            nFrames = dicomFrameCount(info);
            if cvi.contours(c).FrameIndex0 + 1 <= nFrames
                selected = candidates{f};
                selectedInfo = info;
                break;
            end
        end
        if isempty(selected)
            unmatched(c) = true;
            continue;
        end
        frameNumber = cvi.contours(c).FrameIndex0 + 1;
        geometry = cvi_dicom_frame_geometry(selectedInfo,frameNumber);
        points = cvi.contours(c).PointsPixel;
        x = points(:,1).'; % increasing DICOM column index
        y = points(:,2).'; % increasing DICOM row index
        patientLPS = geometry.imagePositionPatient + ...
            geometry.rowDirection .* (geometry.columnSpacingMm .* x) + ...
            geometry.columnDirection .* (geometry.rowSpacingMm .* y);
        % Keep the DICOM patient coordinates in their native LPS-like subject
        % basis.  RS2D IMAGE_ORIENTATION_SUBJECT in these FCI datasets uses
        % the same numerical basis.  Converting to RAH here would place the
        % clinical and FCI geometry in different coordinate systems.
        patientSubject = patientLPS;
        cvi.contours(c).DicomFile = selected;
        cvi.contours(c).PointsPatientMm = patientSubject.';
        cvi.contours(c).FramePositionPatientMm = geometry.imagePositionPatient(:).';
        cvi.contours(c).FrameNormalPatient = geometry.normalDirection(:).';
        cvi.contours(c).SliceIndex = geometry.inStackPositionNumber;
        cvi.contours(c).DicomPhaseIndex = geometry.phaseIndex;
    end

    cvi.unmatchedContourCount = nnz(unmatched);
    cvi.contours = cvi.contours(~unmatched);
    if isempty(cvi.contours)
        error('cvi_import_contours:NoGeometryMatches', ...
            ['Contours were decoded, but none could be matched to image DICOMs. ', ...
             'Select the DICOM export used when the cvi42 workspace was created.']);
    end

    % Series are represented by the enhanced DICOM object containing the
    % relevant frames.  Build phase/surface metadata for each one.
    sourceFiles = unique({cvi.contours.DicomFile},'stable');
    seriesTemplate = struct('DicomFile','','Description','','ContourIndices',[], ...
        'PhaseIndices',[],'PhaseLabels',{{}},'PhaseItemsData',[], ...
        'SourceCenterPatientMm',[0 0 0],'Surfaces',[], ...
        'NumberOfFrames',1,'NumberOfPhases',1,'ContourCount',0, ...
        'ClinicalStackNormalPatient',[NaN NaN NaN],'ClinicalSliceIndices',[]);
    cvi.series = repmat(seriesTemplate,0,1);

    for s = 1:numel(sourceFiles)
        path = sourceFiles{s};
        info = infoCache(path);
        indices = find(strcmp({cvi.contours.DicomFile},path));
        nFrames = dicomFrameCount(info);
        nPhases = cvi_infer_phase_count(info);
        for ii = 1:numel(indices)
            idx = indices(ii);
            % Prefer the enhanced-DICOM frame content.  The supplied Philips
            % cine stores DimensionIndexValues as [stack, slice, cardiacPhase].
            % Fall back to the historical frame-number modulo rule only when
            % the frame metadata are unavailable.
            phase = cvi.contours(idx).DicomPhaseIndex;
            if isempty(phase) || ~isscalar(phase) || ~isfinite(phase) || ...
                    phase < 1 || phase > nPhases
                phase = mod(cvi.contours(idx).FrameIndex0,nPhases) + 1;
            end
            cvi.contours(idx).PhaseIndex = round(phase);
        end
        phases = unique([cvi.contours(indices).PhaseIndex]);
        labels = phaseLabels(cvi.contours,indices,phases);
        series = seriesTemplate;
        series.DicomFile = path;
        if isfield(info,'SeriesDescription') && ~isempty(info.SeriesDescription)
            series.Description = char(string(info.SeriesDescription));
        else
            [~,series.Description] = fileparts(path);
        end
        series.ContourIndices = indices;
        series.PhaseIndices = phases;
        series.PhaseLabels = labels;
        series.PhaseItemsData = phases;
        series.NumberOfFrames = nFrames;
        series.NumberOfPhases = nPhases;
        series.ContourCount = numel(indices);
        series.SourceCenterPatientMm = sourceCentre(cvi.contours,indices);
        series.Surfaces = cvi_build_surfaces(cvi.contours(indices));
        [series.ClinicalStackNormalPatient,series.ClinicalSliceIndices] = ...
            clinicalStackMetadata(cvi.contours(indices));
        cvi.series(end+1,1) = series; %#ok<AGROW>
    end

    cvi.coordinateConvention = 'DICOM/RS2D subject coordinates (LPS-like numerical basis)';

    % Prefer the image object containing the largest number of ventricular
    % contours.  This automatically selects the 330-frame SAX cine in the
    % supplied test study and ignores incidental contours on other series.
    [~,cvi.defaultSeriesIndex] = max([cvi.series.ContourCount]);
end

function count = dicomFrameCount(info)
    count = 1;
    if isfield(info,'NumberOfFrames') && ~isempty(info.NumberOfFrames)
        count = double(info.NumberOfFrames);
        if ischar(info.NumberOfFrames) || isstring(info.NumberOfFrames)
            count = str2double(string(info.NumberOfFrames));
        end
    end
    if ~isfinite(count) || count < 1, count = 1; end
    count = round(count);
end

function centre = sourceCentre(contours,indices)
    preferred = indices(strcmp({contours(indices).Type},'LVEpi'));
    if isempty(preferred)
        preferred = indices(strcmp({contours(indices).Type},'LVEndo'));
    end
    if isempty(preferred), preferred = indices; end
    pointCells = arrayfun(@(x) x.PointsPatientMm,contours(preferred), ...
        'UniformOutput',false);
    points = vertcat(pointCells{:});
    centre = mean(points,1,'omitnan');
end


function [normal,sliceIndices] = clinicalStackMetadata(contours)
    normal = [NaN NaN NaN];
    sliceIndices = [];
    if isempty(contours), return; end
    values = [contours.SliceIndex];
    values = values(isfinite(values));
    if ~isempty(values), sliceIndices = unique(round(values)); end
    normals = vertcat(contours.FrameNormalPatient);
    if ~isempty(normals)
        normals = normals(all(isfinite(normals),2),:);
        if ~isempty(normals)
            candidate = mean(normals,1);
            if norm(candidate) > 0, normal = candidate ./ norm(candidate); end
        end
    end
end

function labels = phaseLabels(contours,indices,phases)
    labels = cell(1,numel(phases));
    meanAreas = nan(1,numel(phases));
    for p = 1:numel(phases)
        mask = strcmp({contours(indices).Type},'LVEndo') & ...
            [contours(indices).PhaseIndex] == phases(p);
        selected = indices(mask);
        if ~isempty(selected)
            meanAreas(p) = mean([contours(selected).PixelArea],'omitnan');
        end
        labels{p} = sprintf('CVI phase %d',phases(p));
    end
    valid = find(isfinite(meanAreas));
    if numel(valid) == 2
        [~,largeLocal] = max(meanAreas(valid));
        [~,smallLocal] = min(meanAreas(valid));
        ed = valid(largeLocal); es = valid(smallLocal);
        labels{ed} = sprintf('ED (CVI phase %d)',phases(ed));
        labels{es} = sprintf('ES (CVI phase %d)',phases(es));
    end
end
