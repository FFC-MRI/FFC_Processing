function [obj, results] = ffc_simple_t1dispersion(obj, slice, options)
%FFC_SIMPLE_T1DISPERSION Simple IR dispersion fit for SVS or an imaging ROI.
%
% Sign selection is performed independently for every evolution field. The
% automatic result for each field is a single early-negative/late-positive
% transition. Candidate transitions are only accepted when the fitted T1 is
% inside a physically plausible range (default 20--1500 ms). User choices in
% the interactive Check-fits window are treated as hard polarity constraints;
% the remaining unconstrained points are re-selected automatically so that a
% plausible T1 is preferred whenever one is available.

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

    plausibleT1RangeMs = [20 1500];
    if isfield(options,'plausibleT1RangeMs') && ...
            isnumeric(options.plausibleT1RangeMs) && ...
            numel(options.plausibleT1RangeMs) == 2
        candidateRange = sort(double(options.plausibleT1RangeMs(:).'));
        if all(isfinite(candidateRange)) && candidateRange(1) > 0 && ...
                candidateRange(2) > candidateRange(1)
            plausibleT1RangeMs = candidateRange;
        end
    end

    fields_mT = double(obj.fieldpoints(:).');
    nFields = numel(fields_mT);
    magnitudeSignals = extractMagnitudeSignals(obj, slice, nFields);
    nTimes = size(magnitudeSignals, 1);
    timesMs = timeMatrix(obj.timepoints, nTimes, nFields);

    % Sign order is defined in increasing evolution time independently for
    % every field, regardless of acquisition/storage ordering.
    sortedTimes = nan(size(timesMs));
    sortedMagnitude = nan(size(magnitudeSignals));
    for fieldIndex = 1:nFields
        [sortedTimes(:, fieldIndex), order] = sort(timesMs(:, fieldIndex));
        sortedMagnitude(:, fieldIndex) = magnitudeSignals(order, fieldIndex);
    end

    fieldTemplate = struct('fieldIndex', NaN, 'field_mT', NaN, ...
        'magnitudeSignal', [], 'automaticSignedSignal', [], ...
        'signedSignal', [], 'timesMs', [], 'prediction', [], ...
        'A', NaN, 'C', NaN, 'm0', NaN, 'T1ms', NaN, 'R1', NaN, ...
        'rawT1ms', NaN, 'rawR1', NaN, ...
        'sse', NaN, 'rsquare', NaN, 'rmse', NaN, 'dfe', NaN, ...
        'automaticNegativePrefixCount', NaN, ...
        'negativePrefixCount', NaN, 'negativeIndices', [], ...
        'positiveIndices', [], 'manualNegativeIndices', [], ...
        'manualPositiveIndices', [], 'manualOverrideCount', 0, ...
        'fitWithinPlausibleRange', false, 'plausibleT1RangeMs', plausibleT1RangeMs, ...
        'candidateT1ms', [], 'candidateScoresUnconstrained', [], ...
        'candidateScoresPlausible', [], 'plausibleCandidateMask', []);
    fieldResults = repmat(fieldTemplate, 1, nFields);
    T1ms = nan(1, nFields);
    R1 = nan(1, nFields);
    selectedNegativePrefixCount = nan(1, nFields);
    candidateScores = inf(max(nTimes - 1, 1), nFields);
    candidateT1ms = nan(max(nTimes - 1, 1), nFields);
    candidatePlausible = false(max(nTimes - 1, 1), nFields);

    for fieldIndex = 1:nFields
        magnitude = sortedMagnitude(:, fieldIndex);
        manualNegative = manualIndicesForField( ...
            obj, 'manual_negative_indices', fieldIndex, nTimes);
        manualPositive = manualIndicesForField( ...
            obj, 'manual_positive_indices', fieldIndex, nTimes);
        overlap = intersect(manualNegative, manualPositive);
        if ~isempty(overlap)
            error('ffc_simple_t1dispersion:ConflictingManualPolarity', ...
                ['Field %d constrains sample(s) %s as both negative and ', ...
                 'positive.'], fieldIndex, mat2str(overlap));
        end

        [prefixCount, fitResult, signed, scoresPlausible, scoresRaw, ...
            t1Candidates, plausibleMask, hasPlausibleFit] = ...
            selectPlausibleSignPrefix(sortedTimes(:, fieldIndex), magnitude, ...
                manualNegative, manualPositive, plausibleT1RangeMs);

        selectedNegativePrefixCount(fieldIndex) = prefixCount;
        candidateScores(1:numel(scoresPlausible),fieldIndex) = scoresPlausible(:);
        candidateT1ms(1:numel(t1Candidates),fieldIndex) = t1Candidates(:);
        candidatePlausible(1:numel(plausibleMask),fieldIndex) = plausibleMask(:);

        automaticSigned = magnitude;
        automaticSigned(1:prefixCount) = -abs(automaticSigned(1:prefixCount));
        automaticSigned((prefixCount + 1):end) = ...
            abs(automaticSigned((prefixCount + 1):end));

        item = fieldTemplate;
        item.fieldIndex = fieldIndex;
        item.field_mT = fields_mT(fieldIndex);
        item.magnitudeSignal = magnitude;
        item.automaticSignedSignal = automaticSigned;
        item.signedSignal = fitResult.signedSignal;
        item.timesMs = fitResult.timesMs;
        item.prediction = fitResult.prediction;
        item.A = fitResult.A;
        item.C = fitResult.C;
        item.m0 = fitResult.m0;
        item.rawT1ms = fitResult.T1ms;
        item.rawR1 = fitResult.R1;
        item.sse = fitResult.sse;
        item.rsquare = fitResult.rsquare;
        item.rmse = fitResult.rmse;
        item.dfe = fitResult.dfe;
        item.fitWithinPlausibleRange = hasPlausibleFit;
        if hasPlausibleFit
            item.T1ms = fitResult.T1ms;
            item.R1 = fitResult.R1;
        else
            % Do not propagate an implausible fit into the dispersion curve.
            % rawT1ms/rawR1 retain the best available fit for diagnostics.
            item.T1ms = NaN;
            item.R1 = NaN;
        end
        item.automaticNegativePrefixCount = prefixCount;
        % Retain the historical field name for compatibility. It now means
        % the selected field-specific prefix after respecting any manual
        % point constraints and the plausible-T1 rule.
        item.negativePrefixCount = prefixCount;
        item.negativeIndices = find(item.signedSignal < 0).';
        item.positiveIndices = find(item.signedSignal >= 0).';
        item.manualNegativeIndices = manualNegative;
        item.manualPositiveIndices = manualPositive;
        item.manualOverrideCount = numel(union(manualNegative,manualPositive));
        item.candidateT1ms = t1Candidates(:).';
        item.candidateScoresUnconstrained = scoresRaw(:).';
        item.candidateScoresPlausible = scoresPlausible(:).';
        item.plausibleCandidateMask = plausibleMask(:).';
        fieldResults(fieldIndex) = item;
        T1ms(fieldIndex) = item.T1ms;
        R1(fieldIndex) = item.R1;
    end

    results = struct();
    results.modelName = 'simple_unweighted_offset_exponential_IR';
    results.modelEquation = 'y(t)=A*exp(-t/T1)+C';
    results.signalDefinition = 'mean magnitude with field-specific signed IR polarity';
    results.fields_mT = fields_mT;
    results.T1ms = T1ms;
    results.R1 = R1;
    results.negativePrefixCount = selectedNegativePrefixCount;
    results.automaticNegativePrefixCount = selectedNegativePrefixCount;
    results.signCandidateScores = candidateScores;
    results.signCandidateT1ms = candidateT1ms;
    results.signCandidatePlausible = candidatePlausible;
    results.plausibleT1RangeMs = plausibleT1RangeMs;
    results.field = fieldResults;
    results.slice = slice;

    obj.dispersioncurve = R1;
    obj.T1FitResults = results;
    obj.param.SimpleIRFit = struct( ...
        'model', results.modelEquation, ...
        'signalDefinition', results.signalDefinition, ...
        'negativePrefixCount', selectedNegativePrefixCount, ...
        'signCandidateScores', candidateScores, ...
        'signCandidateT1ms', candidateT1ms, ...
        'plausibleT1RangeMs', plausibleT1RangeMs, ...
        'magnitudeSignals', sortedMagnitude, ...
        'signedSignals', cat(2, fieldResults.signedSignal));

    if all(selectedNegativePrefixCount == selectedNegativePrefixCount(1))
        prefixText = sprintf('%d', selectedNegativePrefixCount(1));
    else
        prefixText = mat2str(selectedNegativePrefixCount);
    end
    fprintf(['Simple IR fit complete: %d fields, %d time points. ', ...
        'Selected negative prefixes by field: %s. Plausible T1 range %.0f--%.0f ms. ', ...
        '%d/%d fields accepted.\n'], ...
        nFields, nTimes, prefixText, plausibleT1RangeMs(1), ...
        plausibleT1RangeMs(2), sum(isfinite(T1ms)), nFields);

    if optionValue(options, 'makeDispersionPlot', true)
        plotDispersion(obj, results);
    end
    showFits = ~isempty(obj.checkfit) && isscalar(obj.checkfit) && ...
        logical(obj.checkfit);
    if optionValue(options, 'showFits', showFits)
        plotFits(results);
    end
end


function [prefixCount, fitResult, signed, plausibleScores, rawScores, ...
        candidateT1ms, plausibleMask, hasPlausibleFit] = ...
        selectPlausibleSignPrefix(timesMs, magnitude, manualNegative, ...
            manualPositive, plausibleT1RangeMs)
% Select the best single-transition sign pattern, preferring only candidates
% whose fitted T1 lies inside plausibleT1RangeMs. Manual point polarities are
% hard constraints applied to every candidate. If no candidate is plausible,
% retain the lowest-SSE fit for diagnostics but report it as invalid upstream.

    t = double(timesMs(:));
    magnitude = abs(double(magnitude(:)));
    nTimes = numel(magnitude);
    if nTimes < 3
        error('ffc_simple_t1dispersion:TooFewPoints', ...
            'At least three evolution times are required.');
    end

    nCandidates = nTimes - 1;
    rawScores = inf(nCandidates,1);
    plausibleScores = inf(nCandidates,1);
    candidateT1ms = nan(nCandidates,1);
    plausibleMask = false(nCandidates,1);
    fitCache = cell(nCandidates,1);
    signedCache = cell(nCandidates,1);

    scale = sum((magnitude - mean(magnitude,'omitnan')).^2,'omitnan');
    if ~isfinite(scale) || scale <= eps
        scale = max(sum(magnitude.^2,'omitnan'),1);
    end

    for candidate = 1:nCandidates
        trialSigned = magnitude;
        trialSigned(1:candidate) = -abs(trialSigned(1:candidate));
        trialSigned((candidate+1):end) = abs(trialSigned((candidate+1):end));
        if ~isempty(manualNegative)
            trialSigned(manualNegative) = -abs(magnitude(manualNegative));
        end
        if ~isempty(manualPositive)
            trialSigned(manualPositive) = abs(magnitude(manualPositive));
        end
        try
            currentFit = ffc_simple_ir_fit(t,trialSigned);
            fitCache{candidate} = currentFit;
            signedCache{candidate} = trialSigned;
            candidateT1ms(candidate) = currentFit.T1ms;
            rawScores(candidate) = currentFit.sse ./ scale;
            plausibleMask(candidate) = isfinite(currentFit.T1ms) && ...
                currentFit.T1ms >= plausibleT1RangeMs(1) && ...
                currentFit.T1ms <= plausibleT1RangeMs(2);
            if plausibleMask(candidate)
                plausibleScores(candidate) = rawScores(candidate);
            end
        catch
            % Keep this candidate at Inf/NaN and continue searching.
        end
    end

    hasPlausibleFit = any(isfinite(plausibleScores));
    if hasPlausibleFit
        [~,prefixCount] = min(plausibleScores);
    else
        [bestRaw,prefixCount] = min(rawScores);
        if ~isfinite(bestRaw)
            error('ffc_simple_t1dispersion:NoValidSignPattern', ...
                'No inversion-recovery sign pattern could be fitted.');
        end
    end
    fitResult = fitCache{prefixCount};
    signed = signedCache{prefixCount};
end


function indices = manualIndicesForField(obj, propertyName, fieldIndex, nTimes)
    indices = [];
    try
        values = obj.(propertyName);
    catch
        values = [];
    end
    if iscell(values) && fieldIndex <= numel(values) && ...
            ~isempty(values{fieldIndex})
        indices = unique(round(double(values{fieldIndex}(:).')));
        indices = indices(isfinite(indices) & indices >= 1 & indices <= nTimes);
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
    title('Unweighted offset-exponential fit; field-specific IR polarity');
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
        title(ax, sprintf('%.4g mT; T1 %.1f ms; auto -ve %d', ...
            item.field_mT, item.T1ms, item.automaticNegativePrefixCount));
    end
end
