function params = ffc_estimate_tgv_params(images, strength)
%FFC_ESTIMATE_TGV_PARAMS Estimate sensible TGV defaults for an MRI image set.
%
% images:
%   [X Y ...] magnitude image stack.
%
% strength:
%   scalar denoising strength.
%   Suggested:
%       0.5 = light
%       1.0 = moderate
%       2.0 = strong
%       3.0 = very strong
%
% Returns params fields:
%   tgv_alpha0
%   tgv_alpha1
%   tgv_lambda
%   tgv_niter
%   sigma
%   sigmaRel
%   signalScale

if nargin < 2 || isempty(strength)
    strength = 1.0;
end

strength = max(strength, 0.05);

dim = size(images);
stack = abs(reshape(images, dim(1), dim(2), []));

% Use a representative frame
N = size(stack, 3);
k = max(1, round(N/2));
im = single(stack(:,:,k));

% Robust signal scale from non-background-ish pixels
x = double(im(:));
x = x(isfinite(x));

if isempty(x)
    params = fallback_params();
    return
end

% Estimate foreground threshold from robust percentiles
p20 = prctile(x, 20);
p99 = prctile(x, 99);

fg = x > p20 & x < p99;

if nnz(fg) < 100
    fg = x > prctile(x, 50);
end

if nnz(fg) < 100
    fg = x > 0;
end

if nnz(fg) < 100
    signalScale = max(x);
else
    signalScale = median(x(fg));
end

if ~isfinite(signalScale) || signalScale <= 0
    signalScale = max(x);
end

if ~isfinite(signalScale) || signalScale <= 0
    params = fallback_params();
    return
end

% Normalise representative frame
imN = double(im) ./ double(signalScale);

% Estimate noise from high-pass residual using MAD.
% This is not a pure background noise estimate, but is stable for MRI
% magnitude images and does not require a mask.
if exist('imgaussfilt', 'file') == 2
    hp = imN - imgaussfilt(imN, 1.0);
else
    h = fspecial('gaussian', [5 5], 1.0);
    hp = imN - imfilter(imN, h, 'replicate');
end

hpv = hp(:);
hpv = hpv(isfinite(hpv));

sigmaRel = 1.4826 * median(abs(hpv - median(hpv)));

% Clamp to a plausible range so a bad frame does not create absurd params.
sigmaRel = min(max(sigmaRel, 0.002), 0.20);

% Empirical mapping for this TGV2denoise implementation.
% alpha values scale with estimated noise.
params = struct();

params.sigma = sigmaRel * signalScale;
params.sigmaRel = sigmaRel;
params.signalScale = signalScale;

params.tgv_alpha1 = min(max(strength * 2.0 * sigmaRel, 0.005), 0.30);
params.tgv_alpha0 = min(max(strength * 4.0 * sigmaRel, 0.010), 0.60);

% lambda is data fidelity. Larger = less denoising.
% Make lambda decrease as noise/strength increases.
params.tgv_lambda = min(max(1.0 / max(strength * 12.0 * sigmaRel, 0.02), 0.15), 20);

% More iterations for stronger settings.
params.tgv_niter = round(min(max(150 + 80*strength, 100), 500));

end


function params = fallback_params()

params = struct();
params.sigma = [];
params.sigmaRel = [];
params.signalScale = [];

params.tgv_alpha0 = 0.08;
params.tgv_alpha1 = 0.04;
params.tgv_lambda = 1.0;
params.tgv_niter = 250;

end