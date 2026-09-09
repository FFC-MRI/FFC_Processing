function nPhases = cvi_infer_phase_count(info)
%CVI_INFER_PHASE_COUNT Infer cardiac phases in an enhanced cine DICOM.

    candidates = {'CardiacNumberOfImages','NumberOfTemporalPositions'};
    for k = 1:numel(candidates)
        if isfield(info,candidates{k}) && ~isempty(info.(candidates{k}))
            nPhases = scalarNumber(info.(candidates{k}));
            if isfinite(nPhases) && nPhases >= 1
                nPhases = round(nPhases); return;
            end
        end
    end
    nFrames = 1;
    if isfield(info,'NumberOfFrames'), nFrames = scalarNumber(info.NumberOfFrames); end
    nFrames = max(1,round(nFrames));
    if nFrames == 1, nPhases = 1; return; end

    % Fallback: count unique slice positions.  Cine frames from the same
    % slice share position; nFrames / nSlices gives the temporal phase count.
    coords = nan(1,nFrames);
    normal = [];
    for f = 1:nFrames
        try
            g = cvi_dicom_frame_geometry(info,f);
            if isempty(normal), normal = g.normalDirection; end
            coords(f) = dot(normal,g.imagePositionPatient);
        catch
        end
    end
    coords = coords(isfinite(coords));
    if isempty(coords), nPhases = nFrames; return; end
    coords = sort(coords);
    uniqueCoords = coords(1);
    for k = 2:numel(coords)
        if abs(coords(k) - uniqueCoords(end)) > 0.05
            uniqueCoords(end+1) = coords(k); %#ok<AGROW>
        end
    end
    nSlices = numel(uniqueCoords);
    trial = nFrames ./ max(nSlices,1);
    if abs(trial - round(trial)) < 1e-6
        nPhases = max(1,round(trial));
    else
        nPhases = nFrames;
    end
end

function x = scalarNumber(value)
    if ischar(value) || isstring(value), x = str2double(string(value));
    else, x = double(value(1)); end
end
