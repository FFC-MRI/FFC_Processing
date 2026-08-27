function [kspace_selected, info] = prewhiten_and_select(kspace, selectedIndices, whiteningEnabled)
%PREWHITEN_AND_SELECT Prewhiten the full receiver array, then retain a subset.
%
% Receiver data must occupy dimension 6. selectedIndices contains MATLAB's
% one-based indices; the corresponding scanner receiver IDs are therefore
% selectedIndices-1.

    nAcquired = size(kspace, 6);

    if nargin < 2 || isempty(selectedIndices)
        selectedIndices = 1:nAcquired;
    end
    if nargin < 3 || isempty(whiteningEnabled)
        whiteningEnabled = false;
    end

    selectedIndices = unique(round(double(selectedIndices(:).')), 'stable');
    selectedIndices(~isfinite(selectedIndices)) = [];
    if isempty(selectedIndices) || any(selectedIndices < 1) || ...
            any(selectedIndices > nAcquired)
        error('prewhiten_and_select:receiverIndexOutOfRange', ...
            'Selected receiver indices must lie in the range 1..%d.', nAcquired);
    end

    doWhiten = isscalar(whiteningEnabled) && logical(whiteningEnabled) && ...
        nAcquired > 1;
    if doWhiten
        [kspace_prepared, info] = noise_whiten(kspace);
    elseif nAcquired <= 1
        kspace_prepared = kspace;
        info = struct('method', 'singleAcquiredChannel', 'nSamples', 0);
    else
        kspace_prepared = kspace;
        info = struct('method', 'disabled', 'nSamples', 0);
    end

    % Selection occurs only after the full-array whitening transform.
    kspace_selected = kspace_prepared(:,:,:,:,:,selectedIndices);
    info.nInputChannels = nAcquired;
    info.nOutputChannels = numel(selectedIndices);
    info.selectedIndices = selectedIndices;
    info.selectedReceiverIds = selectedIndices - 1;
end
