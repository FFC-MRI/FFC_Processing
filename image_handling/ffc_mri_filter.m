function filtered_images = ffc_mri_filter(images, filter_type, params)
%FFC_MRI_FILTER Denoising pipeline for magnitude MRI image stacks.
%
% images:
%   [X Y ...] magnitude or complex images. Complex input is converted to abs().
%
% filter_type:
%   'None'
%   'Wavelet'
%   'LLR'
%   'LLR+Wavelet'
%   'Total Variation'
%   'Generalised Total Variation'
%   'LLR+Wavelet+TGV'
%
% params:
%   struct with optional fields:
%       patch
%       step
%       niter
%       tau
%       sigma
%       wavelet_level
%       tv_lambda
%       tv_niter
%       tgv_alpha0
%       tgv_alpha1
%       tgv_lambda
%       tgv_niter
%
% Backwards-compatible aliases:
%       tgv_iter    -> tgv_lambda
%       tgv_maxiter -> tgv_niter
%       tv_iter     -> tv_niter
%
% Notes:
%   - All dimensions after X and Y are flattened into the stack dimension.
%   - LLR therefore operates across the flattened image set.
%   - TVL1denoise returns images normalised to [0,1], so this wrapper
%     rescales TV output back to the original image intensity scale.

    if nargin < 2 || isempty(filter_type)
        filter_type = 'None';
    end

    filter_type = normalise_filter_type(filter_type);

    if nargin < 3
        params = [];
    end

    params = normalise_filter_params(params, filter_type);

    dim = size(images);

    if numel(dim) < 2
        error('ffc_mri_filter:InvalidInput', ...
            'Input images must have at least two dimensions.');
    end

    temp = abs(reshape(images, dim(1), dim(2), []));

    inputClass = class(temp);
    temp = single(temp);

    switch lower(filter_type)

        case 'none'
            % passthrough magnitude

        case 'wavelet'
            temp = apply_wavelet(temp, params);

        case 'llr'
            temp = apply_llr(temp, params);

        case 'llr+wavelet'
            temp = apply_llr(temp, params);
            temp = apply_wavelet(temp, params);

        case 'total variation'
            temp = apply_tv(temp, params);

        case 'generalised total variation'
            temp = apply_tgv(temp, params);

        case 'llr+wavelet+tgv'
            temp = apply_llr(temp, params);
            temp = apply_wavelet(temp, params);
            temp = apply_tgv(temp, params);

        otherwise
            warning('ffc_mri_filter:UnknownFilterAfterNormalisation', ...
                'Unknown filter "%s". Returning passthrough magnitude.', filter_type);
    end

    temp = max(temp, 0);

    if strcmp(inputClass, 'double')
        temp = double(temp);
    end

    filtered_images = reshape(temp, dim);

end


% ========================================================================
% Filter type handling
% ========================================================================

function filter_type = normalise_filter_type(filter_type)

    filter_type = char(string(filter_type));
    key = lower(strtrim(filter_type));

    switch key

        case {'none', 'off', 'no filter', 'passthrough', ''}
            filter_type = 'none';

        case {'wavelet', 'wdenoise', 'wavelet denoise'}
            filter_type = 'wavelet';

        case {'llr', 'local low rank', 'locally low rank', 'local low-rank'}
            filter_type = 'llr';

        case {'llr+wavelet', 'llr + wavelet', ...
              'local low rank + wavelet', 'local low-rank + wavelet'}
            filter_type = 'llr+wavelet';

        case {'total variation', 'tv', 'tvl1', 'tv-l1', 'tv l1'}
            filter_type = 'total variation';

        case {'generalised total variation', 'generalized total variation', ...
              'tgv', 'tgv2', 'tgv2denoise'}
            filter_type = 'generalised total variation';

        case {'llr+wavelet+tgv', 'llr + wavelet + tgv', ...
              'llr+wavelet+generalised total variation', ...
              'llr + wavelet + generalised total variation'}
            filter_type = 'llr+wavelet+tgv';

        case {'deep learning'}
            % Historical GUI option. No implementation here.
            filter_type = 'none';

        otherwise
            warning('ffc_mri_filter:UnknownFilter', ...
                'Unknown filter type "%s". Using passthrough.', filter_type);
            filter_type = 'none';
    end

end


% ========================================================================
% Parameter handling
% ========================================================================

function params = normalise_filter_params(params, filter_type)

    defaults = default_filter_params(filter_type);

    if isempty(params) || ~isstruct(params)
        params = defaults;
        return
    end

    % Backwards-compatible aliases from earlier GUI versions.
    if isfield(params, 'tgv_iter') && ~isfield(params, 'tgv_lambda')
        params.tgv_lambda = params.tgv_iter;
    end

    if isfield(params, 'tgv_maxiter') && ~isfield(params, 'tgv_niter')
        params.tgv_niter = params.tgv_maxiter;
    end

    if isfield(params, 'tv_iter') && ~isfield(params, 'tv_niter')
        params.tv_niter = params.tv_iter;
    end

    % If the GUI only has one lambda/iteration control, allow the TGV
    % lambda/iteration values to drive TV too.
    if isfield(params, 'tgv_lambda') && ~isfield(params, 'tv_lambda')
        params.tv_lambda = params.tgv_lambda;
    end

    if isfield(params, 'tgv_niter') && ~isfield(params, 'tv_niter')
        params.tv_niter = params.tgv_niter;
    end

    fn = fieldnames(defaults);

    for k = 1:numel(fn)
        f = fn{k};

        if ~isfield(params, f) || isempty(params.(f))
            params.(f) = defaults.(f);
            continue
        end

        if isnumeric(params.(f)) && any(~isfinite(double(params.(f)(:))))
            params.(f) = defaults.(f);
        end
    end

    params.patch        = round(params.patch);
    params.step         = round(params.step);
    params.niter        = round(params.niter);
    params.wavelet_level = round(params.wavelet_level);
    params.tv_niter     = round(params.tv_niter);
    params.tgv_niter    = round(params.tgv_niter);

    params.patch = max(params.patch, 3);
    if mod(params.patch, 2) == 0
        params.patch = params.patch + 1;
    end

    params.step = max(params.step, 1);
    params.niter = max(params.niter, 1);
    params.wavelet_level = max(params.wavelet_level, 1);

    params.tau = max(params.tau, 0);

    params.tv_lambda = max(params.tv_lambda, eps);
    params.tv_niter = max(params.tv_niter, 1);

    params.tgv_alpha0 = max(params.tgv_alpha0, 0);
    params.tgv_alpha1 = max(params.tgv_alpha1, 0);
    params.tgv_lambda = max(params.tgv_lambda, eps);
    params.tgv_niter = max(params.tgv_niter, 1);

end


function params = default_filter_params(filter_type)

    params = struct();

    % LLR defaults
    params.patch = 5;
    params.step = 3;
    params.niter = 2;
    params.tau = 0.5;
    params.sigma = [];

    % Wavelet defaults
    params.wavelet_level = 3;

    % TV defaults
    % TVL1denoise returns [0,1], and lambda interpretation can be quite
    % strong. Start conservative.
    params.tv_lambda = 0.20;
    params.tv_niter = 100;

    % TGV defaults
    params.tgv_alpha0 = 0.04;
    params.tgv_alpha1 = 0.015;
    params.tgv_lambda = 25;
    params.tgv_niter = 150;

    switch lower(filter_type)

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

        case 'wavelet'
            params.wavelet_level = 3;

        case 'total variation'
            params.tv_lambda = 0.20;
            params.tv_niter = 100;

        case 'generalised total variation'
            params.tgv_alpha0 = 0.04;
            params.tgv_alpha1 = 0.015;
            params.tgv_lambda = 25;
            params.tgv_niter = 150;

        case 'llr+wavelet+tgv'
            params.patch = 8;
            params.step = 3;
            params.niter = 2;
            params.tau = 1.2;
            params.wavelet_level = 3;
            params.tgv_alpha0 = 0.04;
            params.tgv_alpha1 = 0.015;
            params.tgv_lambda = 25;
            params.tgv_niter = 150;
    end

end


% ========================================================================
% Individual filter wrappers
% ========================================================================

function temp = apply_wavelet(temp, params)

    if exist('wavelet_per_frame', 'file') ~= 2
        warning('ffc_mri_filter:MissingWaveletFunction', ...
            'wavelet_per_frame.m was not found. Returning input unchanged.');
        return
    end

    temp = wavelet_per_frame(temp, params.wavelet_level);

end


function temp = apply_llr(temp, params)

    if exist('llr_denoise_stack', 'file') ~= 2
        warning('ffc_mri_filter:MissingLLRFunction', ...
            'llr_denoise_stack.m was not found. Returning input unchanged.');
        return
    end

    temp = llr_denoise_stack(temp, struct( ...
        'patch', params.patch, ...
        'step', params.step, ...
        'niter', params.niter, ...
        'tau', params.tau, ...
        'sigma', params.sigma));

end


function temp = apply_tv(temp, params)

    if exist('TVL1denoise', 'file') ~= 2
        warning('ffc_mri_filter:MissingTVL1Function', ...
            'TVL1denoise.m was not found. Falling back to TGV.');
        temp = apply_tgv(temp, params);
        return
    end

    tmp = zeros(size(temp), 'like', temp);

    if size(temp, 3) == 1
        tmp(:,:,1) = tvl1_preserve_scale(temp(:,:,1), ...
            params.tv_lambda, params.tv_niter);
    else
        parfor n = 1:size(temp, 3)
            tmp(:,:,n) = tvl1_preserve_scale(temp(:,:,n), ...
                params.tv_lambda, params.tv_niter);
        end
    end

    temp = tmp;

end


function out = tvl1_preserve_scale(im, lambda, niter)

    im = single(im);

    mx = max(im(:));

    if ~isfinite(mx) || mx <= 0
        out = im;
        return
    end

    % TVL1denoise normalises internally and returns an image in [0,1].
    % Rescale to the original dynamic range so GUI windowing still works.
    out = TVL1denoise(im, lambda, niter);

    out = single(out) .* mx;
    out = max(out, 0);

end


function temp = apply_tgv(temp, params)

    if exist('TGV2denoise', 'file') ~= 2
        warning('ffc_mri_filter:MissingTGVFunction', ...
            'TGV2denoise.m was not found. Returning input unchanged.');
        return
    end

    optsT = struct('data', 'L1', 'normalize', true);

    tmp = zeros(size(temp), 'like', temp);

    if size(temp, 3) == 1
        tmp(:,:,1) = TGV2denoise(temp(:,:,1), ...
            params.tgv_alpha0, ...
            params.tgv_alpha1, ...
            params.tgv_lambda, ...
            params.tgv_niter, ...
            optsT);
    else
        parfor n = 1:size(temp, 3)
            tmp(:,:,n) = TGV2denoise(temp(:,:,n), ...
                params.tgv_alpha0, ...
                params.tgv_alpha1, ...
                params.tgv_lambda, ...
                params.tgv_niter, ...
                optsT);
        end
    end

    temp = tmp;

end