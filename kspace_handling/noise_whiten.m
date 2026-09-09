function [out, info] = noise_whiten(kspace, opts)
%NOISE_WHITEN Robust calibration plus bounded full-channel whitening.
% Public API retains channels in the last dimension, as before.
    if nargin < 2, opts = struct(); end
    info = ffc_noise_calibration(kspace, opts);
    sz = size(kspace);
    if numel(sz) < 3 || sz(end) <= 1, out = kspace; return; end
    nCoils = sz(end);
    input = reshape(kspace, sz(1), sz(2), [], 1, 1, nCoils);
    out = ffc_apply_whitening(input, info, 1:nCoils, nCoils);
    out = reshape(out, sz);
end
