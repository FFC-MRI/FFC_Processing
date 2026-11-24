function [plan, skipped] = ffc_plan_dicom_batch(scans, folder)
% Snapshot each scan's own metadata on the client; never share selected patient.
    plan = struct('index',{},'filename',{},'options',{},'patientKey',{});
    skipped = {};
    metadata = cell(size(scans)); errors = cell(size(scans));
    for k = 1:numel(scans)
        try, metadata{k} = ffc_extract_dicom_metadata(scans{k});
        catch ME, errors{k} = ME; end
    end
    for k = 1:numel(scans)
        scan = scans{k};
        try
            if isempty(scan.magimage) || (~isempty(scan.is_svs) && logical(scan.is_svs))
                skipped{end+1} = sprintf('%d: no imaging magnitude data (SVS excluded)',k); %#ok<AGROW>
                continue;
            end
            if ~isempty(errors{k}), rethrow(errors{k}); end
            meta = metadata{k};
            % Fill only missing demographics, with unanimous values from the
            % same explicit PatientID. Never borrow from the selected scan or
            % match anonymous/missing IDs by position in the loaded list.
            if ~isempty(meta.PatientID)
                for field = {'PatientName','PatientBirthDate','PatientSex'}
                    name = field{1};
                    if ~isempty(meta.(name)), continue; end
                    values = {};
                    for q = 1:numel(metadata)
                        other = metadata{q};
                        if ~isempty(other) && strcmp(meta.PatientID,other.PatientID) && ~isempty(other.(name))
                            values{end+1} = other.(name); %#ok<AGROW>
                        end
                    end
                    values = unique(values);
                    if numel(values) == 1, meta.(name) = values{1}; end
                end
            end
            if isempty(meta.PatientID) && isempty(meta.PatientName)
                error('ffc:MissingIdentity','No patient ID or name: cannot safely assign this series to a participant.');
            end
            if ~isempty(meta.PatientID)
                patientKey = jsonencode({'ID',meta.PatientID});
            else
                patientKey = jsonencode({'NAME',meta.PatientName,meta.PatientBirthDate});
            end
            % Keep separate participants AND separate visits/studies distinct.
            studyKey = jsonencode({patientKey,meta.StudyID,meta.StudyDescription,meta.StudyDate});
            if isempty(meta.StudyDate)
                % Without a visit date, avoid assuming unrelated scans share a study.
                studyKey = jsonencode({studyKey,meta.AcquisitionUUID,k});
            end
            number = str2double(meta.SeriesID);
            if ~isscalar(number) || ~isfinite(number) || number <= 0, number = 9000+k; end
            label = meta.PatientID; if isempty(label), label = meta.PatientName; end
            filename = sprintf('%05d_%s_%s_EnhancedMR.dcm',k,safeName(label,40), ...
                safeName(meta.SeriesDescription,60));
            opts = struct('ResolvedMetadata',meta, ...
                'StudyInstanceUID',ffc_make_dicom_uid(['FCI-BATCH-STUDY|' studyKey]), ...
                'FrameOfReferenceUID',ffc_make_dicom_uid(['FCI-BATCH-FOR-LPS25|' studyKey]), ...
                'SeriesInstanceUID',ffc_make_dicom_uid(), ...
                'SOPInstanceUID',ffc_make_dicom_uid(),'SeriesNumber',round(number));
            plan(end+1) = struct('index',k,'filename',fullfile(folder,filename), ...
                'options',opts,'patientKey',patientKey); %#ok<AGROW>
        catch ME
            skipped{end+1} = sprintf('%d: %s',k,ME.message); %#ok<AGROW>
        end
    end
end
function value = safeName(value,maxLength)
    value = regexprep(char(string(value)),'[^A-Za-z0-9_-]+','_');
    if isempty(value), value = 'Series'; end
    value = value(1:min(maxLength,numel(value)));
end
