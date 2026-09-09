function result = ffc_compute_task(kind, input, payload)
%FFC_COMPUTE_TASK Pure computation shared by background and batch execution.
    switch kind
        case 'preview'
            payload.params.use_parallel = false;
            result = ffc_mri_filter(input,payload.filter,payload.params);
        case 'reconstruct'
            originalParams = input.denoise_params;
            if isempty(input.denoise_params), input.denoise_params = struct(); end
            input.denoise_params.use_parallel = false;
            result = buildimages(input,payload.force);
            result.denoise_params = originalParams;
        case 'map'
            result = ffc_simple_t1map(input,struct('use_parallel',false));
        case 'filter'
            result = input;
            source = input.filter_source_magimage;
            if isempty(source), source = input.magimage; end
            params = payload.params; params.use_parallel = false;
            result.magimage = ffc_mri_filter(source,payload.filter,params);
            result.filter_source_magimage = source;
            if ~isempty(result.is_svs) && result.is_svs
                result.param.SVS_PROCESSED_MAGNITUDE_SIGNAL = reshape(result.magimage, ...
                    size(result.magimage,4),size(result.magimage,5));
            end
            result.denoise_filter = payload.filter;
            result.denoise_params = payload.params;
            result.registration_tform = [];
            result = ffc_invalidate_analysis(result);
            result.processing_history{end+1} = struct('operation','filter', ...
                'filter',payload.filter,'params',payload.params);
        case 'register'
            result = input;
            result.magimage = ffc_register_images(input.magimage,payload.masks, ...
                struct('use_parallel',false));
            % A later filter starts from this registered result. Rebuild
            % returns to native reconstruction and clears registration.
            result.filter_source_magimage = result.magimage;
            result.denoise_filter = 'None';
            result.denoise_params = ffc_default_denoise_params('None');
            result = ffc_invalidate_analysis(result);
            result.registration_tform = struct('method','demons', ...
                'masks',payload.masks,'iterations',[400 200 100]);
            result.processing_history{end+1} = struct('operation','register', ...
                'iterations',[400 200 100]);
        otherwise
            error('ffc_compute_task:UnknownTask','Unknown task: %s',kind);
    end
end
