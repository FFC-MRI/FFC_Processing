function [scan, timing] = ffc_open_cameleon_scan(folder, openRoot, reprocess)
%FFC_OPEN_CAMELEON_SCAN Reuse baseline images before reading raw data.
    if nargin < 3, reprocess = false; end
    timing = struct('readSeconds',0,'cacheLookupSeconds',0,'cacheHit',false);
    cacheTimer = tic;
    scan = [];
    if ~reprocess, scan = ffc_disk_cache('load',folder); end
    timing.cacheLookupSeconds = toc(cacheTimer);
    if ~isempty(scan)
        timing.cacheHit = true;
        scan.source_metadata.OpenRootPath = openRoot;
        return;
    end
    readTimer = tic;
    before = ffc_disk_cache('identity',folder);
    [data,par] = read_raw_data_Cameleon(folder,[]);
    par.openRootPath = openRoot;
    scan = ImageReconCore(struct('data',data,'par',par),'.dat',1);
    scan.disk_cache_source = folder;
    scan.disk_cache_identity = before;
    timing.readSeconds = toc(readTimer);
end
