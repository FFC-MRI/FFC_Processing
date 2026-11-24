function report = ffc_export_dicom_batch(scans, folder, progress)
% Export every loaded imaging scan with bounded parallelism, in a fresh folder.
% Progress callback returns false to cancel. Completed files are retained.
    if nargin < 3, progress = []; end
    started = tic;
    token = strrep(ffc_make_dicom_uid(),'.','_');
    runFolder = fullfile(folder,['FCI_export_' token]);
    [ok,message] = mkdir(runFolder);
    if ~ok, error('ffc:ExportFolder','%s',message); end
    [plan,skipped] = ffc_plan_dicom_batch(scans,runFolder);
    report = struct('folder',runFolder,'results',{{}},'skipped',{skipped}, ...
        'failed',{{}},'cancelled',false,'workers',0,'seconds',0,'participantCount',0);
    if ~isempty(plan), report.participantCount = numel(unique({plan.patientKey})); end
    count = numel(plan); completed = 0; next = 1;
    jobs = {}; jobIndices = []; pool = [];
    futureScope = FFCFutureScope();
    cleanup = futureScope.guard(); %#ok<NASGU>
    limit = getpref('FFCProcessing','ExportWorkers',2);
    validateattributes(limit,{'numeric'},{'scalar','integer','positive','finite'});
    if count > 1 && getpref('FFCProcessing','ParallelDicomExport',true)
        if ~update('Starting export workers'), report.cancelled=true; return; end
        try
            if ~license('test','Distrib_Computing_Toolbox')
                error('ffc:NoParallelToolbox','Parallel Computing Toolbox unavailable.');
            end
            pool = gcp('nocreate');
            if isempty(pool), pool = parpool('local',min(limit,count)); end
            if isa(pool,'parallel.ThreadPool')
                error('ffc:ThreadPool','Existing thread pool retained; using serial export.');
            end
            limit = min(limit,pool.NumWorkers); report.workers = limit;
        catch ME
            pool = [];
            warning('ffc:SerialExport','Using serial DICOM export: %s',ME.message);
        end
    end
    while next <= count || ~isempty(jobs)
        % Collect finished jobs before testing cancellation, keeping their reports.
        for j = numel(jobs):-1:1
            if ~strcmp(jobs{j}.State,'finished'), continue; end
            try
                result = fetchOutputs(jobs{j});
                report.results{end+1} = result;
            catch ME
                report.failed{end+1} = sprintf('%d: %s',plan(jobIndices(j)).index,ME.message);
            end
            completed = completed+1;
            futureScope.release(jobs{j});
            delete(jobs{j}); jobs(j)=[]; jobIndices(j)=[];
        end
        if ~update(sprintf('Exported/processed %d of %d series',completed,count))
            report.cancelled = true; break;
        end
        if next <= count && numel(jobs) < limit
            n = next; next = next+1;
            if isempty(pool)
                try
                    report.results{end+1} = ffc_write_dicom_job(exportInput(scans{plan(n).index}),plan(n));
                catch ME
                    report.failed{end+1} = sprintf('%d: %s',plan(n).index,ME.message);
                end
                completed = completed+1;
            else
                try
                    job = parfeval(pool,@ffc_write_dicom_job,1,exportInput(scans{plan(n).index}),plan(n));
                    futureScope.add(job);
                    jobs{end+1}=job; jobIndices(end+1)=n; %#ok<AGROW>
                catch ME
                    report.failed{end+1} = sprintf('%d: %s',plan(n).index,ME.message);
                    completed = completed+1;
                end
            end
        elseif ~isempty(jobs)
            pause(0.05);
        end
    end
    futureScope.cancelAll();
    report.seconds = toc(started);
    if ~report.cancelled, update('Export complete'); end
    function proceed = update(label)
        proceed = true;
        if ~isempty(progress), proceed = progress(label,completed,count); end
    end
end

function scan = exportInput(scan)
% For valid scans only, transfer magnitude + metadata rather than raw/coil cubes.
% This is a value-class snapshot; the viewer retains its full reconstruction.
    if ~ffc_reconstruction_valid(scan), return; end
    for name = {'rawdata','originalcomplexkspace','complexkspace','compleximage', ...
            'complexcombined','phaseimage','magkspace','filter_source_magimage', ...
            'sensitivity_cache','whitening_cache','T1Maps','R1Maps'}
        scan.(name{1}) = [];
    end
end
