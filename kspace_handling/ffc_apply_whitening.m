function out = ffc_apply_whitening(in, info, selected, nCoils, blockSamples)
%FFC_APPLY_WHITENING Bounded scratch space, full acquired-channel transform.
% Selection is AFTER whitening, including when only one output is requested.
    if nargin < 5, blockSamples = 65536; end
    sz = size(in); sz(end+1:6) = 1;
    nSamples = numel(in) / nCoils;
    input = reshape(in, nSamples, nCoils);
    if ~isfield(info, 'psi') || isempty(info.psi)
        output = input(:, selected);
    else
        output = zeros(nSamples, numel(selected), 'like', in);
        for first = 1:blockSamples:nSamples
            rows = first:min(nSamples, first + blockSamples - 1);
            X = input(rows, :).';
            if info.p == 0
                X = info.L \ X;
            else
                X = info.Winv * X;
            end
            output(rows, :) = X(selected, :).';
        end
    end
    % The reconstruction API explicitly keeps receivers in dimension six.
    sz(6) = numel(selected);
    out = reshape(output, sz);
end
