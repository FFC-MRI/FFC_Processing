function [obj, results] = ffc_simple_t1dispersion(obj, slice, options)
%FFC_SIMPLE_T1DISPERSION Simple IR dispersion fit for SVS or an imaging ROI.

    if nargin < 2 || isempty(slice)
        slice = 1;
    end
    if nargin < 3 || isempty(options)
        options = struct();
    end
    if ~isstruct(options)
        error('ffc_simple_t1dispersion:InvalidOptions', ...
            'Options must be a struct.');
    end

    fields_mT = double(obj.fieldpoints(:).');
    nFields = numel(fields_mT);
    magnitudeSignals = extractMagnitudeSignals(obj, slice, nFields);
    nTimes = size(magnitudeSignals, 1);
    timesMs = timeMatrix(obj.timepoints, nTimes, nFields);

    % Sign order is defined in increasing evolution time, independently of
    % the order in which fields or experiments were stored.
    sortedTimes = nan(size(timesMs));
    sortedMagnitude = nan(size(magnitudeSignals));
    for fieldIndex = 1:nFields
        [sortedTimes(:, fieldIndex), order] = ...
            sort(timesMs(:, fieldIndex));
        sortedMagnitude(:, fieldIndex) = ...
            magnitudeSignals(order, fieldIndex);
    end

    [negativePrefixCount, candidateScores, selectedFits] = ...
        ffc_select_ir_prefix(sortedTimes, sortedMagnitude);

    fieldTemplate = struct('fieldIndex', NaN, 'field_mT', NaN, ...
        'magnitudeSignal', [], 'signedSignal', [], 'timesMs', [], ...
        'prediction', [], 'A', NaN, 'C', NaN, 'm0', NaN, ...
        'T1ms', NaN, 'R1', NaN, 'sse', NaN, 'rsquare', NaN, ...
        'rmse', NaN, 'dfe', NaN, 'negativePrefixCount', NaN);
    fieldResults = repmat(fieldTemplate, 1, nFields);
    T1ms = nan(1, nFields);
    R1 = nan(1, nFields);
    for fieldIndex = 1:nFields
        fitResult = selectedFits{fieldIndex};
        item = fieldTemplate;
        item.fieldIndex = fieldIndex;
        item.field_mT = fields_mT(fieldIndex);
        item.magnitudeSignal = sortedMagnitude(:, fieldIndex);
        item.signedSignal = fitResult.signedSignal;
        item.timesMs = fitResult.timesMs;
        item.prediction = fitResult.prediction;
        item.A = fitResult.A;
        item.C = fitResult.C;
        item.m0 = fitResult.m0;
        item.T1ms = fitResult.T1ms;
        item.R1 = fitResult.R1;
        item.sse = fitResult.sse;
        item.rsquare = fitResult.rsquare;
        item.rmse = fitResult.rmse;
        item.dfe = fitResult.dfe;
        item.negativePrefixCount = negativePrefixCount;
        fieldResults(fieldIndex) = item;
        T1ms(fieldIndex) = item.T1ms;
        R1(fieldIndex) = item.R1;
    end

    results = struct();
    results.modelName = 'simple_unweighted_offset_exponential_IR';
    results.modelEquation = 'y(t)=A*exp(-t/T1)+C';
    results.signalDefinition = 'mean magnitude';
    results.fields_mT = fields_mT;
    results.T1ms = T1ms;
    results.R1 = R1;
    results.negativePrefixCount = negativePrefixCount;
    results.signCandidateScores = candidateScores;
    results.field = fieldResults;
    results.slice = slice;

    obj.dispersioncurve = R1;
    obj.T1FitResults = results;
    obj.param.SimpleIRFit = struct( ...
        'model', results.modelEquation, ...
        'signalDefinition', results.signalDefinition, ...
        'negativePrefixCount', negativePrefixCount, ...
        'signCandidateScores', candidateScores, ...
        'magnitudeSignals', sortedMagnitude, ...
        'signedSignals', cat(2, fieldResults.signedSignal));

    fprintf(['Simple IR fit complete: %d fields, %d time points, ', ...
        '%d data-selected early negative points.\n'], ...
        nFields, nTimes, negativePrefixCount);

    if optionValue(options, 'makeDispersionPlot', true)
        plotDispersion(obj, results);
    end
    showFits = ~isempty(obj.checkfit) && isscalar(obj.checkfit) && ...
        logical(obj.checkfit);
    if optionValue(options, 'showFits', showFits)
        plotFits(results);
    end
end


function signals = extractMagnitudeSignals(obj, slice, nFields)
    isSVS = ~isempty(obj.is_svs) && isscalar(obj.is_svs) && ...
        logical(obj.is_svs);
    if isSVS && isstruct(obj.param) && ...
            isfield(obj.param, 'SVS_PROCESSED_MAGNITUDE_SIGNAL')
        signals = double(obj.param.SVS_PROCESSED_MAGNITUDE_SIGNAL);
        if size(signals, 2) ~= nFields && size(signals, 1) == nFields
            signals = signals.';
        end
        if size(signals, 2) ~= nFields
            error('ffc_simple_t1dispersion:SVSFieldMismatch', ...
                'The SVS signal does not contain all %d fields.', nFields);
        end
        return;
    end

    data = obj.magimage;
    dataSize = size(data);
    dataSize(end + 1:5) = 1;
    if slice < 1 || slice > dataSize(3)
        error('ffc_simple_t1dispersion:SliceOutOfRange', ...
            'Slice %d is outside 1..%d.', slice, dataSize(3));
    end
    if dataSize(5) ~= nFields
        error('ffc_simple_t1dispersion:FieldMismatch', ...
            'Magnitude data contain %d fields; fieldpoints contains %d.', ...
            dataSize(5), nFields);
    end
    mask = obj.mask;
    if isempty(mask)
        mask = true(dataSize(1), dataSize(2));
    elseif ndims(mask) >= 3
        mask = mask(:, :, min(slice, size(mask, 3)));
    end
    mask = logical(mask);
    if ~isequal(size(mask), dataSize(1:2)) || ~any(mask(:))
        error('ffc_simple_t1dispersion:InvalidMask', ...
            'The selected imaging ROI is empty or has the wrong size.');
    end

    nTimes = dataSize(4);
    signals = nan(nTimes, nFields);
    for fieldIndex = 1:nFields
        for timeIndex = 1:nTimes
            frame = data(:, :, slice, timeIndex, fieldIndex);
            values = double(frame(mask & isfinite(frame)));
            if ~isempty(values)
                signals(timeIndex, fieldIndex) = mean(abs(values));
            end
        end
    end
end


function matrix = timeMatrix(timepoints, nTimes, nFields)
    values = double(timepoints);
    if isvector(values)
        if numel(values) == nTimes
            matrix = repmat(values(:), 1, nFields);
        elseif numel(values) == nTimes .* nFields
            matrix = reshape(values, [nTimes, nFields]);
        else
            error('ffc_simple_t1dispersion:TimeCountMismatch', ...
                'Cannot map %d time values to %d x %d data.', ...
                numel(values), nTimes, nFields);
        end
    elseif isequal(size(values), [nFields, nTimes])
        matrix = values.';
    elseif isequal(size(values), [nTimes, nFields])
        matrix = values;
    else
        error('ffc_simple_t1dispersion:TimeSizeMismatch', ...
            'timepoints has size %s; expected [%d %d] or [%d %d].', ...
            mat2str(size(values)), nFields, nTimes, nTimes, nFields);
    end
end


function value = optionValue(options, name, defaultValue)
    value = defaultValue;
    if isfield(options, name) && isscalar(options.(name))
        value = logical(options.(name));
    end
end


function plotDispersion(obj, results)
    plotR1 = ~isempty(obj.R1T1) && isscalar(obj.R1T1) && logical(obj.R1T1);
    figure('Name', 'Simple inversion-recovery dispersion');
    if plotR1
        plot(results.fields_mT, results.R1, 'o-', 'LineWidth', 1.2);
        ylabel('R_1 (s^{-1})');
    else
        plot(results.fields_mT, results.T1ms, 'o-', 'LineWidth', 1.2);
        ylabel('T_1 (ms)');
    end
    xlabel('Evolution field (mT)');
    title(sprintf(['Unweighted offset-exponential fit; ', ...
        '%d early points selected negative'], ...
        results.negativePrefixCount));
    grid on;
end


function plotFits(results)
    nFields = numel(results.field);
    nColumns = ceil(sqrt(nFields));
    nRows = ceil(nFields ./ nColumns);
    figure('Name', 'Simple inversion-recovery fits');
    layout = tiledlayout(nRows, nColumns, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for fieldIndex = 1:nFields
        item = results.field(fieldIndex);
        ax = nexttile(layout);
        plot(ax, item.timesMs, item.signedSignal, 'o');
        hold(ax, 'on');
        denseTime = linspace(min(item.timesMs), max(item.timesMs), 200).';
        plot(ax, denseTime, item.A .* exp(-denseTime ./ item.T1ms) + ...
            item.C, '-', 'LineWidth', 1.1);
        yline(ax, 0, '--');
        hold(ax, 'off');
        grid(ax, 'on');
        title(ax, sprintf('%.4g mT; T1 %.1f ms', ...
            item.field_mT, item.T1ms));
    end
end
