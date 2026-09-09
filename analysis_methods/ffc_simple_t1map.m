function obj = ffc_simple_t1map(obj, options)
%FFC_SIMPLE_T1MAP Apply the simple unweighted IR fit to every selected voxel.

    if nargin < 2, options = struct(); end
    magnitudeData = obj.magimage;
    dataSize = size(magnitudeData);
    dataSize(end + 1:5) = 1;
    nTimes = dataSize(4);
    nFields = dataSize(5);
    if numel(obj.fieldpoints) ~= nFields
        error('ffc_simple_t1map:FieldMismatch', ...
            'Magnitude data and fieldpoints contain different field counts.');
    end
    timesMs = mapTimeMatrix(obj.timepoints, nTimes, nFields);

    mask = mapMask(obj, magnitudeData, dataSize);
    maskIndices = find(mask(:));
    if isempty(maskIndices)
        error('ffc_simple_t1map:EmptyMask', ...
            'No voxels were selected for relaxation mapping.');
    end

    flattened = reshape(magnitudeData, ...
        [prod(dataSize(1:3)), nTimes, nFields]);
    globalSignals = nan(nTimes, nFields);
    for fieldIndex = 1:nFields
        for timeIndex = 1:nTimes
            values = double(flattened(maskIndices, timeIndex, fieldIndex));
            globalSignals(timeIndex, fieldIndex) = ...
                mean(abs(values), 'omitnan');
        end
    end

    sortedTimes = nan(size(timesMs));
    sortedSignals = nan(size(globalSignals));
    sortOrder = zeros(nTimes, nFields);
    for fieldIndex = 1:nFields
        [sortedTimes(:, fieldIndex), sortOrder(:, fieldIndex)] = ...
            sort(timesMs(:, fieldIndex));
        sortedSignals(:, fieldIndex) = ...
            globalSignals(sortOrder(:, fieldIndex), fieldIndex);
    end
    [negativePrefixCount, candidateScores] = ...
        ffc_select_ir_prefix(sortedTimes, sortedSignals);

    nSelected = numel(maskIndices);
    selectedT1 = nan(nSelected, nFields);
    selectedR1 = nan(nSelected, nFields);
    selectedRsquare = nan(nSelected, nFields);
    useParallel = ffc_use_parallel(options,nSelected);
    % Pack only a bounded block of selected voxels. parfor sees a sliced input.
    for first = 1:4096:nSelected
        rows = first:min(nSelected,first+4095);
        packed = double(flattened(maskIndices(rows),:,:));
        blockT1 = nan(numel(rows),nFields);
        blockR1 = blockT1; blockRsq = blockT1;
        if useParallel
            parfor i = 1:numel(rows)
                [blockT1(i,:),blockR1(i,:),blockRsq(i,:)] = fitMapVoxel( ...
                    packed(i,:,:),sortedTimes,sortOrder,negativePrefixCount,nFields);
            end
        else
            for i = 1:numel(rows)
                [blockT1(i,:),blockR1(i,:),blockRsq(i,:)] = fitMapVoxel( ...
                    packed(i,:,:),sortedTimes,sortOrder,negativePrefixCount,nFields);
            end
        end
        selectedT1(rows,:) = blockT1; selectedR1(rows,:) = blockR1;
        selectedRsquare(rows,:) = blockRsq;
    end

    fullT1 = nan(prod(dataSize(1:3)), nFields);
    fullR1 = nan(prod(dataSize(1:3)), nFields);
    fullRsquare = nan(prod(dataSize(1:3)), nFields);
    fullT1(maskIndices, :) = selectedT1;
    fullR1(maskIndices, :) = selectedR1;
    fullRsquare(maskIndices, :) = selectedRsquare;
    obj.T1Maps = reshape(fullT1, [dataSize(1:3), nFields]);
    obj.R1Maps = reshape(fullR1, [dataSize(1:3), nFields]);
    obj.param.SimpleIRMapFit = struct( ...
        'model', 'y(t)=A*exp(-t/T1)+C', ...
        'signalDefinition', 'magnitude at each voxel', ...
        'negativePrefixCount', negativePrefixCount, ...
        'signCandidateScores', candidateScores, ...
        'Rsquare', reshape(fullRsquare, [dataSize(1:3), nFields]), ...
        'mask', mask);
    fprintf(['Simple IR mapping complete: %d voxels x %d fields; ', ...
        '%d data-selected early negative points.\n'], ...
        nSelected, nFields, negativePrefixCount);
end


function [T1row, R1row, rsquareRow] = fitMapVoxel( ...
        voxelData, sortedTimes, sortOrder, negativePrefixCount, nFields)
    T1row = nan(1, nFields);
    R1row = nan(1, nFields);
    rsquareRow = nan(1, nFields);
    voxelData = reshape(double(voxelData), [], nFields);
    for fieldIndex = 1:nFields
        magnitude = abs(voxelData(sortOrder(:, fieldIndex), fieldIndex));
        signed = magnitude;
        signed(1:negativePrefixCount) = ...
            -signed(1:negativePrefixCount);
        try
            [T1row(fieldIndex),R1row(fieldIndex),rsquareRow(fieldIndex)] = ...
                ffc_simple_ir_fit(sortedTimes(:,fieldIndex),signed);
        catch
        end
    end
end


function mask = mapMask(obj, magnitudeData, dataSize)
    mask = [];
    if ~isempty(obj.mask)
        candidate = logical(obj.mask);
        if isequal(size(candidate), dataSize(1:3))
            mask = candidate;
        elseif isequal(size(candidate), dataSize(1:2))
            mask = repmat(candidate, [1 1 dataSize(3)]);
        end
    end
    if isempty(mask)
        reference = max(magnitudeData, [], 4);
        reference = max(reference, [], 5);
        threshold = 0.02 .* max(reference(:));
        mask = isfinite(reference) & reference > threshold;
    end
end


function matrix = mapTimeMatrix(timepoints, nTimes, nFields)
    values = double(timepoints);
    if isvector(values)
        if numel(values) == nTimes
            matrix = repmat(values(:), 1, nFields);
        elseif numel(values) == nTimes .* nFields
            matrix = reshape(values, [nTimes, nFields]);
        else
            error('ffc_simple_t1map:TimeCountMismatch', ...
                'The number of evolution times does not match the images.');
        end
    elseif isequal(size(values), [nFields, nTimes])
        matrix = values.';
    elseif isequal(size(values), [nTimes, nFields])
        matrix = values;
    else
        error('ffc_simple_t1map:TimeSizeMismatch', ...
            'The evolution-time matrix does not match the image dimensions.');
    end
end
