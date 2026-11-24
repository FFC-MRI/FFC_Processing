function result = ffc_write_dicom_job(scan, job)
% One independent series per worker; atomically publish a complete DICOM.
    ffc_setup();
    rebuilt = ~ffc_reconstruction_valid(scan);
    if rebuilt
        scan = ffc_compute_task('reconstruct',scan,struct('force',false));
    end
    folder = fileparts(job.filename);
    temporary = [tempname(folder) '.part'];
    cleanup = onCleanup(@() removePart(temporary)); %#ok<NASGU>
    info = ffc_export_enhanced_mr(scan,temporary,job.options);
    if isfile(job.filename)
        error('ffc:ExportCollision','Output already exists: %s',job.filename);
    end
    [ok,message] = movefile(temporary,job.filename);
    if ~ok, error('ffc:ExportMove','%s',message); end
    info.filename = job.filename;
    result = struct('index',job.index,'info',info);
    % Avoid transferring all image arrays back when no rebuild was required.
    result.rebuiltScan = [];
    if rebuilt, result.rebuiltScan = scan; end
end
function removePart(file)
    if isfile(file), delete(file); end
end
