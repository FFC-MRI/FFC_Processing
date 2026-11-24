function [scans, names, failures] = ffc_preprocess_scans(inputScans, options)
%FFC_PREPROCESS_SCANS Prepare a batch without retaining partially loaded scans.
% Only successfully preprocessed scans and their matching display names are
% returned. Each failure retains the original input index and MException so
% the caller can report the full stack. No GUI or existing dataset is changed.

    if nargin < 2
        options = struct();
    end
    scans = cell(1, numel(inputScans));
    names = cell(1, numel(inputScans));
    failures = struct('sequenceIndex', {}, 'exception', {});
    nAccepted = 0;

    for sequenceIndex = 1:numel(inputScans)
        if isfield(options,'progress') && ~isempty(options.progress)
            if ~options.progress(sequenceIndex,numel(inputScans))
                error('ffc:LoadCancelled','Dataset loading cancelled.');
            end
        end
        try
            scan = rename_sequence(inputScans{sequenceIndex});
            if isfield(options, 'backgroundselect')
                scan.backgroundselect = options.backgroundselect;
            end
            if isfield(options, 'reprocess')
                scan.reprocess = options.reprocess;
            end
            if isfield(options, 'checkfit')
                scan.checkfit = options.checkfit;
            end
            scan = preprocessing(scan);
            if isfield(options,'reconstruct') && options.reconstruct
                scan = ffc_normalize_viewer_scan(scan);
                reused = ffc_reconstruction_valid(scan);
                scan = buildimages(scan,false);
                % Only raw-folder baseline reconstructions are disk cached.
                % Manual postprocessing remains in the current session.
                if ~reused && ~isempty(scan.disk_cache_source)
                    ffc_disk_cache('save',scan.disk_cache_source,scan);
                end
            end

            nAccepted = nAccepted + 1;
            scans{nAccepted} = scan;
            names{nAccepted} = scan.displayname;
        catch ME
            failures(end + 1) = struct( ...
                'sequenceIndex', sequenceIndex, 'exception', ME); %#ok<AGROW>
        end
    end

    scans = scans(1:nAccepted);
    names = names(1:nAccepted);
end
