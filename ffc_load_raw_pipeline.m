function [scans, report, failures] = ffc_load_raw_pipeline(folders, openRoot, options)
% One client reader, at most two in-flight reconstruction jobs by default.
% Results are kept in source order. Workers never open the source/cache files.
% Cancellation throws; caller retains the previously displayed dataset.
    if nargin < 3, options = struct(); end
    reader = @ffc_open_cameleon_scan;
    if isfield(options,'reader'), reader = options.reader; end % client-side test/reader seam
    count = numel(folders); started = tic;
    scans = cell(1,count);
    failures = struct('sequenceIndex',{},'exception',{});
    row = struct('readSeconds',0,'cacheLookupSeconds',0,'preprocessSeconds',0, ...
        'reconstructSeconds',0,'cacheWriteSeconds',0,'cacheHit',false,'status','pending');
    report = struct('sequence',repmat(row,1,count),'wallSeconds',0, ...
        'poolStartupSeconds',0,'maxInFlight',0,'workers',0,'fallbackReason','');
    limit = getpref('FFCProcessing','LoadWorkers',2);
    validateattributes(limit,{'numeric'},{'scalar','integer','positive','finite'});
    limit = min(limit,max(1,count));
    parallelEnabled = getpref('FFCProcessing','ParallelLoading',true) && count > 1;
    pool = []; poolChecked = false;
    jobs = {}; indices = []; next = 1; completed = 0;
    futureScope = FFCFutureScope();
    cleanup = futureScope.guard(); %#ok<NASGU>
    workerOptions = struct();
    for name = {'backgroundselect','reprocess','checkfit'}
        if isfield(options,name{1}), workerOptions.(name{1}) = options.(name{1}); end
    end
    reprocess = isfield(options,'reprocess') && options.reprocess;
    while next <= count || ~isempty(jobs)
        progress('Reading / reconstructing');
        % Collect finished work first, freeing queue slots and worker outputs.
        for j = numel(jobs):-1:1
            if ~strcmp(jobs{j}.State,'finished'), continue; end
            index = indices(j);
            try
                result = fetchOutputs(jobs{j});
                accept(result,index);
            catch ME
                failed(index,ME);
            end
            futureScope.release(jobs{j});
            delete(jobs{j}); jobs(j) = []; indices(j) = [];
        end
        if next <= count && numel(jobs) < limit
            index = next; next = next+1;
            progress(sprintf('Reading sequence %d of %d',index,count));
            try
                [scan,timing] = reader(folders{index},openRoot,reprocess);
                report.sequence(index).readSeconds = timing.readSeconds;
                report.sequence(index).cacheLookupSeconds = timing.cacheLookupSeconds;
                report.sequence(index).cacheHit = timing.cacheHit;
            catch ME
                failed(index,ME); continue;
            end
            progress(sprintf('Preparing sequence %d of %d',index,count));
            % Warm cache must not pay pool startup / worker transfer overhead.
            if ~reprocess && ffc_reconstruction_valid(ffc_normalize_viewer_scan(scan))
                try
                    accept(ffc_prepare_loaded_scan(scan,workerOptions),index);
                catch ME
                    failed(index,ME);
                end
                clear scan
                continue;
            end
            if ~poolChecked
                poolChecked = true; poolTimer = tic;
                if parallelEnabled
                    progress('Starting reconstruction workers');
                    try
                        if ~license('test','Distrib_Computing_Toolbox')
                            error('ffc:NoParallelToolbox','Parallel Computing Toolbox unavailable.');
                        end
                        pool = gcp('nocreate');
                        if isempty(pool), pool = parpool('local',limit); end
                        % Process workers isolate toolbox state and support legacy routines.
                        if isa(pool,'parallel.ThreadPool')
                            error('ffc:ThreadPool','Existing thread pool retained; using serial loading.');
                        end
                        limit = min(limit,pool.NumWorkers);
                        report.workers = limit;
                    catch ME
                        pool = []; report.fallbackReason = ME.message;
                        warning('ffc:SerialLoading','Using serial reconstruction: %s',ME.message);
                    end
                end
                report.poolStartupSeconds = toc(poolTimer);
            end
            if isempty(pool)
                try
                    accept(ffc_prepare_loaded_scan(scan,workerOptions),index);
                catch ME
                    failed(index,ME);
                end
            else
                try
                    job = parfeval(pool,@ffc_prepare_loaded_scan,1,scan,workerOptions);
                    futureScope.add(job);
                    jobs{end+1} = job; indices(end+1) = index; %#ok<AGROW>
                    report.maxInFlight = max(report.maxInFlight,numel(jobs));
                catch ME
                    % Submission failure: preserve functionality on the client.
                    report.fallbackReason = ME.message;
                    try
                        accept(ffc_prepare_loaded_scan(scan,workerOptions),index);
                    catch serialError
                        failed(index,serialError);
                    end
                end
            end
            clear scan
        elseif ~isempty(jobs)
            pause(0.05); % yields UI events; no blocking fetchNext
        end
    end
    progress('All sequences processed');
    report.wallSeconds = toc(started);
    report.successCount = nnz(~cellfun(@isempty,scans));
    report.failureCount = numel(failures);
    scans = scans(~cellfun(@isempty,scans));
    if ~isempty(failures)
        [~,order] = sort([failures.sequenceIndex]); failures = failures(order);
    end
    fprintf(['FFC opening: %.1f s elapsed; %.1f s raw/header reads; %.1f s cache lookup; ' ...
        '%.1f s preprocessing; %.1f s reconstruction; %.1f s cache writes; %.1f s pool startup; %d workers.\n'], ...
        report.wallSeconds,sum([report.sequence.readSeconds]), ...
        sum([report.sequence.cacheLookupSeconds]),sum([report.sequence.preprocessSeconds]), ...
        sum([report.sequence.reconstructSeconds]),sum([report.sequence.cacheWriteSeconds]),report.poolStartupSeconds,report.workers);
    % Summed worker times overlap and need not equal elapsed wall time.
    function progress(label)
        if isfield(options,'loadProgress') && ~isempty(options.loadProgress)
            state = struct('label',label,'completed',completed,'total',count, ...
                'active',numel(jobs),'elapsed',toc(started), ...
                'readSeconds',sum([report.sequence.readSeconds]), ...
                'reconstructSeconds',sum([report.sequence.reconstructSeconds]));
            if ~options.loadProgress(state)
                error('ffc:LoadCancelled','Dataset loading cancelled.');
            end
        end
    end
    function accept(result,index)
        report.sequence(index).preprocessSeconds = result.preprocessSeconds;
        report.sequence(index).reconstructSeconds = result.reconstructSeconds;
        cacheTimer = tic;
        if ~result.reused && ~isempty(result.scan.disk_cache_source)
            ffc_disk_cache('save',result.scan.disk_cache_source,result.scan);
        end
        report.sequence(index).cacheWriteSeconds = toc(cacheTimer);
        scans{index} = result.scan;
        report.sequence(index).status = 'ready';
        completed = completed+1;
    end
    function failed(index,exception)
        failures(end+1) = struct('sequenceIndex',index,'exception',exception);
        report.sequence(index).status = 'failed'; completed = completed+1;
    end
end
