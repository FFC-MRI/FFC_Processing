function [complexSignal, perReceiverSignal, magnitudeSignal, ...
        perReceiverMagnitudeSignal] = svs_signal_from_kspace(data)
%SVS_SIGNAL_FROM_KSPACE Reduce a localised acquisition to one signal curve.
%
% DATA follows the pipeline convention
%   [sample, phase, partition, time, field, receiver].
% SVS must have singleton phase and partition dimensions.  The returned
% COMPLEXSIGNAL and MAGNITUDESIGNAL are [1 1 1 time field]. The former is
% retained for phase display/QC. The quantitative relaxometry signal follows
% the validated direct extraction: MEAN(ABS(DATA),1), so magnitude is taken
% before averaging sampled points. Multiple selected receivers are averaged
% after their per-receiver reductions; this changes only overall scale.

    if isempty(data)
        error('svs_signal_from_kspace:EmptyData', ...
            'SVS data cannot be empty.');
    end
    if size(data, 2) ~= 1 || size(data, 3) ~= 1
        error('svs_signal_from_kspace:NotSingleVoxel', ...
            ['Expected singleton phase and partition dimensions, but got ', ...
             '[%d %d].'], size(data, 2), size(data, 3));
    end

    nTimes = size(data, 4);
    nFields = size(data, 5);
    nReceivers = size(data, 6);

    perReceiverSignal = mean(data, 1);
    perReceiverSignal = reshape(perReceiverSignal, ...
        [1, 1, 1, nTimes, nFields, nReceivers]);
    complexSignal = mean(perReceiverSignal, 6);
    complexSignal = reshape(complexSignal, ...
        [1, 1, 1, nTimes, nFields]);

    perReceiverMagnitudeSignal = mean(abs(data), 1);
    perReceiverMagnitudeSignal = reshape(perReceiverMagnitudeSignal, ...
        [1, 1, 1, nTimes, nFields, nReceivers]);
    magnitudeSignal = mean(perReceiverMagnitudeSignal, 6);
    magnitudeSignal = reshape(magnitudeSignal, ...
        [1, 1, 1, nTimes, nFields]);
end
