function filtered_images = ffc_mri_filter(images, filter_type, params)
%FFC_MRI_FILTER Denoising pipeline for magnitude stacks.
%
% images      : [X Y ...] magnitude or complex images; magnitude is used.
% filter_type : 'None', 'Wavelet', 'LLR', 'LLR+Wavelet',
%               'Total Variation', 'Generalised Total Variation',
%               'LLR+Wavelet+TGV'
% params      : struct with optional fields
%               patch, step, niter, tau, sigma, wavelet_level,
%               tgv_alpha0, tgv_alpha1, tgv_lambda, tgv_niter

if nargin < 2 || isempty(filter_type)
    filter_type = 'None';
end

if nargin < 3 || isempty(params) || ~isstruct(params)
    params = ffc_default_denoise_params(filter_type);
else
    defaults = ffc_default_denoise_params(filter_type);
    fn = fieldnames(defaults);
    for k = 1:numel(fn)
        if ~isfield(params, fn{k}) || isempty(params.(fn{k}))
            params.(fn{k}) = defaults.(fn{k});
        end
    end
end

if mod(params.patch,2) == 0
    params.patch = params.patch + 1;
end

filter_key = lower(strtrim(filter_type));

dim = size(images);
temp = abs(reshape(images, dim(1), dim(2), []));

switch filter_key
    case {'none','off','deep learning'}
        % passthrough

    case 'llr+wavelet'
        temp = llr_denoise_stack(temp, struct( ...
            'patch', params.patch, 'step', params.step, ...
            'niter', params.niter, 'tau', params.tau, 'sigma', params.sigma));
        temp = wavelet_per_frame(temp, params.wavelet_level);

    case 'llr'
        temp = llr_denoise_stack(temp, struct( ...
            'patch', params.patch, 'step', params.step, ...
            'niter', params.niter, 'tau', params.tau, 'sigma', params.sigma));

    case 'wavelet'
        temp = wavelet_per_frame(temp, params.wavelet_level);

    case {'total variation','generalised total variation','llr+wavelet+tgv'}
        temp = llr_denoise_stack(temp, struct( ...
            'patch', params.patch, 'step', params.step, ...
            'niter', params.niter, 'tau', params.tau, 'sigma', params.sigma));
        temp = wavelet_per_frame(temp, params.wavelet_level);

        optsT = struct('data','L1','normalize',true);
        tmp2 = zeros(size(temp), 'like', temp);
        parfor n = 1:size(temp,3)
            tmp2(:,:,n) = TGV2denoise(temp(:,:,n), ...
                params.tgv_alpha0, params.tgv_alpha1, ...
                params.tgv_lambda, params.tgv_niter, optsT);
        end
        temp = tmp2;

    otherwise
        % passthrough
end

filtered_images = reshape(temp, dim);
end
