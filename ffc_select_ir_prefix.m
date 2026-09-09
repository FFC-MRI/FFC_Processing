function [negativePrefixCount, candidateScores, selectedFits] = ...
        ffc_select_ir_prefix(timesMs, magnitudeSignals)
%FFC_SELECT_IR_PREFIX Select one data-driven IR sign transition.
%
% Every candidate has at least one early negative sample and one late
% positive sample. The same candidate index is assessed across all fields,
% matching an acquisition in which every field uses the same ordered IR
% schedule. Each curve is fitted by ordinary unweighted least squares and
% contributes a scale-normalised SSE to the selection score.

    timesMs = double(timesMs);
    magnitudeSignals = double(magnitudeSignals);
    if ~isequal(size(timesMs), size(magnitudeSignals))
        error('ffc_select_ir_prefix:SizeMismatch', ...
            'Times and magnitude signals must have the same size.');
    end
    [nTimes, nFields] = size(magnitudeSignals);
    if nTimes < 3
        error('ffc_select_ir_prefix:TooFewPoints', ...
            'At least three time points are required for an IR fit.');
    end

    candidateScores = inf(1, nTimes - 1);
    fitCache = cell(nTimes - 1, nFields);
    for candidate = 1:(nTimes - 1)
        totalScore = 0;
        for fieldIndex = 1:nFields
            magnitude = magnitudeSignals(:, fieldIndex);
            signed = magnitude;
            signed(1:candidate) = -signed(1:candidate);
            try
                currentFit = ffc_simple_ir_fit( ...
                    timesMs(:, fieldIndex), signed);
                fitCache{candidate, fieldIndex} = currentFit;
                scale = sum((magnitude - mean(magnitude, 'omitnan')).^2, ...
                    'omitnan');
                if ~isfinite(scale) || scale <= eps
                    scale = max(sum(magnitude.^2, 'omitnan'), 1);
                end
                totalScore = totalScore + currentFit.sse ./ scale;
            catch
                totalScore = Inf;
                break;
            end
        end
        candidateScores(candidate) = totalScore;
    end

    [bestScore, negativePrefixCount] = min(candidateScores);
    if ~isfinite(bestScore)
        error('ffc_select_ir_prefix:NoValidCandidate', ...
            'No single-transition inversion-recovery fit succeeded.');
    end
    selectedFits = fitCache(negativePrefixCount, :);
end
