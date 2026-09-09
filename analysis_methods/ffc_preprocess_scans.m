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
