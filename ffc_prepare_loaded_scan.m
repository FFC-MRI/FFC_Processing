function result = ffc_prepare_loaded_scan(scan, options)
% Pure worker stage: no source reads, cache writes or UI access.
    ffc_setup();
    scan = rename_sequence(scan);
    for name = {'backgroundselect','reprocess','checkfit'}
        if isfield(options,name{1}), scan.(name{1}) = options.(name{1}); end
    end
    started = tic;
    scan = preprocessing(scan);
    preSeconds = toc(started);
    scan = ffc_normalize_viewer_scan(scan);
    reused = ffc_reconstruction_valid(scan);
    started = tic;
    scan = ffc_compute_task('reconstruct',scan,struct('force',false));
    result = struct('scan',scan,'preprocessSeconds',preSeconds, ...
        'reconstructSeconds',toc(started),'reused',reused);
end
