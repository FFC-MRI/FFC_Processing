function params = ffc_default_denoise_params(filter_type)
%FFC_DEFAULT_DENOISE_PARAMS Default parameter struct for the denoiser GUI.

if nargin < 1 || isempty(filter_type)
    filter_type = 'LLR+Wavelet';
end

filter_type = char(string(filter_type));
filter_key = lower(strtrim(filter_type));

params = struct();

% LLR defaults
params.patch = 5;
params.step = 3;
params.niter = 2;
params.tau = 0.5;
params.sigma = [];

% Wavelet defaults
params.wavelet_level = 3;

% Total variation defaults
params.tv_lambda = 0.20;
params.tv_niter = 100;

% TGV defaults
params.tgv_alpha0 = 0.04;
params.tgv_alpha1 = 0.015;
params.tgv_lambda = 25;
params.tgv_niter = 150;

switch filter_key

    case {'none','off','deep learning'}
        % passthrough; keep generic defaults

    case 'wavelet'
        params.wavelet_level = 3;

    case 'llr'
        params.patch = 8;
        params.step = 3;
        params.niter = 2;
        params.tau = 2.2;

    case 'llr+wavelet'
        params.patch = 5;
        params.step = 3;
        params.niter = 2;
        params.tau = 0.5;
        params.wavelet_level = 3;

    case {'total variation','tv','tvl1','tv-l1'}
        params.tv_lambda = 0.20;
        params.tv_niter = 100;

    case {'generalised total variation','generalized total variation','tgv'}
        params.tgv_alpha0 = 0.04;
        params.tgv_alpha1 = 0.015;
        params.tgv_lambda = 25;
        params.tgv_niter = 150;

    case {'llr+wavelet+tgv','llr + wavelet + tgv'}
        params.patch = 8;
        params.step = 3;
        params.niter = 2;
        params.tau = 1.2;
        params.wavelet_level = 3;

        params.tgv_alpha0 = 0.04;
        params.tgv_alpha1 = 0.015;
        params.tgv_lambda = 25;
        params.tgv_niter = 150;

    otherwise
        % Keep generic defaults.
end

end