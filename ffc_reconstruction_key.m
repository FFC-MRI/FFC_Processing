function key = ffc_reconstruction_key(scan)
%FFC_RECONSTRUCTION_KEY Small, explicit signature; never hashes image arrays.
% External scripts replacing raw/preprocessed data must call preprocessing
% with reprocess=true, or buildimages(scan) (the default is a forced rebuild).
    names = {'source_generation','samples','views','slices','n_timepoints', ...
        'n_fieldpoints','n_receivers','fft_size','fft_size_3d', ...
        'window_size','window_function','multichannel_recon','combination', ...
        'biascorrect','noise_whitening','TwoDimensional','is_svs', ...
        'orientation','fov','fov_phase','fov_3d','thk','retain_coil_images','retain_phase_images'};
    key = struct();
    key.algorithmVersion = 31;
    key.rfCorrection = ffc_rf_correction_config(scan);
    key.rf_phase_order_override = scan.rf_phase_order_override;
    for i = 1:numel(names)
        key.(names{i}) = scan.(names{i});
    end
    % Only parameters actually consumed by numerical reconstruction belong
    % here. Analysis, export and contour metadata never invalidate images.
    parameters = {'OFF_CENTER_FIELD_OF_VIEW_1D','OFF_CENTER_FIELD_OF_VIEW_2D', ...
        'OFF_CENTER_FIELD_OF_VIEW_3D','OFF_CENTER_FIELD_OF_VIEW_EFF', ...
        'FIELD_OF_VIEW','FIELD_OF_VIEW_PHASE','FIELD_OF_VIEW_3D'};
    key.param = struct();
    if isstruct(scan.param)
        for i = 1:numel(parameters)
            name = parameters{i};
            if isfield(scan.param,name), key.param.(name) = scan.param.(name); end
        end
    elseif iscell(scan.param)
        % Legacy parameter rows use the same names as correct_orientation.
        selected = false(size(scan.param,1),1);
        for i = 1:size(scan.param,1)
            selected(i) = any(strcmpi(strtrim(string(scan.param{i,1})),parameters));
        end
        key.param = scan.param(selected,:);
    end
end
