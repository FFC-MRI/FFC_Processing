function cvi = cvi_import_contours(workspaceFile,dicomRoot)
%CVI_IMPORT_CONTOURS Import cvi42 contours and attach clinical DICOM geometry.
%
%   CVI = CVI_IMPORT_CONTOURS(WORKSPACEFILE,DICOMROOT) reads ventricular
%   contours and converts every matched point into DICOM LPS millimetres.
%   TEST25 FCI geometry is also returned in LPS by ffc_scan_spatial_geometry.
%   Clinical DICOM points must not receive the scanner Y/Z conversion.

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

    decodedContours = cvi.contours;
    directMatched = false(numel(decodedContours),1);
    for c = 1:numel(decodedContours)
        uidIndex = find(strcmp(uniqueUIDs,decodedContours(c).ImageUID),1);
        candidates = uidMatches{uidIndex};
        selected = '';
        selectedInfo = [];
        for f = 1:numel(candidates)
            info = infoCache(candidates{f});
            if isempty(info) || ~isfield(info,'Rows') || ~isfield(info,'Columns')
                continue;
            end
            nFrames = dicomFrameCount(info);
            if decodedContours(c).FrameIndex0 + 1 <= nFrames
                selected = candidates{f};
                selectedInfo = info;
                break;
            end
        end
        if isempty(selected), continue; end
        frameNumber = decodedContours(c).FrameIndex0 + 1;
        geometry = cvi_dicom_frame_geometry(selectedInfo,frameNumber);
        points = decodedContours(c).PointsPixel;
        x = points(:,1).'; % increasing DICOM column index
        y = points(:,2).'; % increasing DICOM row index
        patientLPS = geometry.imagePositionPatient + ...
            geometry.rowDirection .* (geometry.columnSpacingMm .* x) + ...
            geometry.columnDirection .* (geometry.rowSpacingMm .* y);
        decodedContours(c).DicomFile = selected;
        decodedContours(c).PointsPatientMm = patientLPS.';
        decodedContours(c).FramePositionPatientMm = geometry.imagePositionPatient(:).';
        decodedContours(c).FrameNormalPatient = geometry.normalDirection(:).';
        decodedContours(c).SliceIndex = geometry.inStackPositionNumber;
        decodedContours(c).DicomPhaseIndex = geometry.phaseIndex;
        directMatched(c) = true;
    end

    % cvi42 can create derived processing images (notably T1 maps) with its
    % own 1.3.6.1.4.1.53684... UIDs, then copy the source contours onto
    % those images. Scanner exports need not contain those derived images.
    % Detect exact contour-set duplicates and treat them as aliases instead
    % of reporting them as missing geometry or inserting duplicate surfaces.
    [derivedDuplicate,duplicateGroups] = detectDerivedDuplicateContours(decodedContours,directMatched);
    unresolved = ~directMatched & ~derivedDuplicate;
    primaryMask = ismember({decodedContours.Type},{'LVEndo','LVEpi','RVEndo'}).';
    cvi.decodedContourCount = numel(decodedContours);
    cvi.directMatchedContourCount = nnz(directMatched);
    cvi.derivedDuplicateContourCount = nnz(derivedDuplicate);
    cvi.derivedDuplicateVentricularCount = nnz(derivedDuplicate & primaryMask);
    cvi.unmatchedContourCount = nnz(unresolved);
    cvi.unmatchedVentricularContourCount = nnz(unresolved & primaryMask);
    cvi.derivedDuplicateGroups = duplicateGroups;
    cvi.contours = decodedContours(directMatched);
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
        'ClinicalStackNormalPatient',[NaN NaN NaN],'ClinicalSliceIndices',[], ...
        'PrimaryContourCount',0,'LGEEnhancementAvailable',false, ...
        'LGEOverlayCount',0,'LGEOverlayDetails',struct());
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
        primary = ismember({cvi.contours(indices).Type},{'LVEndo','LVEpi','RVEndo'});
        series.PrimaryContourCount = nnz(primary);
        series.SourceCenterPatientMm = sourceCentre(cvi.contours,indices);
        series.Surfaces = cvi_build_surfaces(cvi.contours(indices));

        % The cvi42 LGE coloured enhancement overlay is not a normal contour
        % object. It is stored as a sparse WorkImageOverlay mask. Reconstruct
        % it and attach a patient-space mesh to the matching LGE series.
        [lgeSurface,lgeDetails] = cvi_build_lge_enhancement_surface( ...
            cvi.workImageOverlays,cvi.contours(indices),info);
        series.LGEEnhancementAvailable = lgeDetails.available;
        series.LGEOverlayCount = lgeDetails.matchedOverlayCount;
        series.LGEOverlayDetails = lgeDetails;
        if lgeDetails.available
            series.Surfaces(end+1,1) = lgeSurface; %#ok<AGROW>
        end
        [series.ClinicalStackNormalPatient,series.ClinicalSliceIndices] = ...
            clinicalStackMetadata(cvi.contours(indices));
        cvi.series(end+1,1) = series; %#ok<AGROW>
    end

    cvi.coordinateConvention = 'DICOM LPS (TEST25 shared patient geometry)';

    % When a saved cvi42 LGE enhancement overlay is present, select that
    % series by default because it is the most information-rich contour set.
    % Otherwise preserve the original behaviour of choosing the series with
    % the largest ventricular contour stack.
    lge = find([cvi.series.LGEEnhancementAvailable]);
    if ~isempty(lge)
        [~,q] = max([cvi.series(lge).LGEOverlayCount]);
        cvi.defaultSeriesIndex = lge(q);
    else
        [~,cvi.defaultSeriesIndex] = max([cvi.series.PrimaryContourCount]);
    end
end

function [isDuplicate,groups] = detectDerivedDuplicateContours(contours,directMatched)
%DETECTDERIVEDDUPLICATECONTOURS Identify cvi42-only derived image aliases.
% Require an exact LV endocardial + epicardial pair match. If an ancillary
% blood-pool contour is also present it must match as well. This is much
% safer than aliasing a derived UID from a single circular ROI alone.

    isDuplicate = false(numel(contours),1);
    groups = struct('DerivedImageKey',{},'SourceImageKey',{},'ContourCount',{});
    allKeys = {contours.ImageKey};
    unmatchedKeys = unique(allKeys(~directMatched),'stable');
    matchedKeys = unique(allKeys(directMatched),'stable');
    for u = 1:numel(unmatchedKeys)
        key = unmatchedKeys{u};
        idx = find(strcmp(allKeys,key));
        if isempty(idx), continue; end
        typesHere = {contours(idx).Type};
        if ~all(ismember({'LVEndo','LVEpi'},typesHere)), continue; end
        candidates = {};
        for m = 1:numel(matchedKeys)
            sourceKey = matchedKeys{m};
            src = find(strcmp(allKeys,sourceKey) & directMatched.');
            sourceTypes = {contours(src).Type};
            required = {'LVEndo','LVEpi'};
            if any(strcmp(typesHere,'BloodPool')), required{end+1} = 'BloodPool'; end %#ok<AGROW>
            ok = true;
            for t = 1:numel(required)
                a = idx(strcmp(typesHere,required{t}));
                b = src(strcmp(sourceTypes,required{t}));
                if numel(a) ~= 1 || numel(b) ~= 1 || ...
                        ~samePoints(contours(a).PointsPixel,contours(b).PointsPixel)
                    ok = false; break;
                end
            end
            if ok, candidates{end+1} = sourceKey; end %#ok<AGROW>
        end
        if numel(candidates) ~= 1, continue; end
        sourceKey = candidates{1};
        src = find(strcmp(allKeys,sourceKey) & directMatched.');
        sourceTypes = {contours(src).Type};
        % Mark only contour types that have an exact counterpart on the
        % source image. Unknown derived-only objects remain unresolved.
        marked = false(size(idx));
        for q = 1:numel(idx)
            b = src(strcmp(sourceTypes,contours(idx(q)).Type));
            if numel(b)==1 && samePoints(contours(idx(q)).PointsPixel,contours(b).PointsPixel)
                isDuplicate(idx(q)) = true; marked(q) = true;
            end
        end
        if any(marked)
            g.DerivedImageKey = key;
            g.SourceImageKey = sourceKey;
            g.ContourCount = nnz(marked);
            groups(end+1,1) = g; %#ok<AGROW>
        end
    end
end

function tf = samePoints(a,b)
    tf = isequal(size(a),size(b));
    if ~tf, return; end
    if isempty(a), tf = true; return; end
    tf = all(abs(double(a(:))-double(b(:))) < 1e-9);
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
