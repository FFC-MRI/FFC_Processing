function key = ffc_reconstruction_key(scan)
%FFC_RECONSTRUCTION_KEY Small, explicit signature; never hashes image arrays.
% External scripts replacing raw/preprocessed data must call preprocessing
% with reprocess=true, or buildimages(scan) (the default is a forced rebuild).
    names = {'source_generation','samples','views','slices','n_timepoints', ...
        'n_fieldpoints','n_receivers','fft_size','fft_size_3d', ...
        'window_size','window_function','multichannel_recon','combination', ...
        'biascorrect','noise_whitening','TwoDimensional','is_svs', ...
        'orientation','fov','fov_phase','fov_3d','thk','param','retain_coil_images'};
    key = struct();
    for i = 1:numel(names)
        key.(names{i}) = scan.(names{i});
    end
    % SVS fitting adds diagnostics to param; those do not change the source.
    if isstruct(key.param)
        fields = fieldnames(key.param);
        diagnostics = startsWith(fields, 'SVS_PROCESSED_') | ...
            startsWith(fields, 'SVS_PER_RECEIVER_') | startsWith(fields,'SimpleIR');
        key.param = rmfield(key.param, fields(diagnostics));
    end
end
