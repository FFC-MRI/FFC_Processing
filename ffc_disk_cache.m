function output = ffc_disk_cache(action, folder, scan)
%FFC_DISK_CACHE Local baseline reconstruction cache, never writes source data.
% Raw identity uses size/mtime; header and Serie.xml are compared byte-for-byte.
% Force Reprocessing to bypass after source edits preserving raw size/mtime.
    output = [];
    if ~getpref('FFCProcessing','DiskCacheEnabled',true), return; end
    % Allow reuse of existing caches without blocking loading on new writes.
    if strcmp(action,'save') && ~getpref('FFCProcessing','DiskCacheWriteEnabled',true), return; end
    try
        jf = java.io.File(folder);
        folder = char(jf.getCanonicalPath());
        stamp = identity(folder);
        if isempty(stamp), return; end
        if strcmp(action,'identity'), output = stamp; return; end
        cacheRoot = getpref('FFCProcessing','DiskCacheDirectory', ...
            fullfile(prefdir,'FFC_Reconstruction_Cache_TEST26'));
        digest = java.security.MessageDigest.getInstance('SHA-256');
        digest.update(unicode2native(folder,'UTF-8'));
        token = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
        target = fullfile(cacheRoot,[token '.mat']);
        if strcmp(action,'load')
            if ~isfile(target), return; end
            meta = load(target,'cacheIdentity');
            if ~isfield(meta,'cacheIdentity') || ~isequaln(meta.cacheIdentity,stamp), return; end
            data = load(target,'cachedScan');
            if isa(data.cachedScan,'ImageReconCore') && ffc_reconstruction_valid(data.cachedScan) && ...
                    isequaln(stamp,identity(folder))
                output = data.cachedScan;
                output.disk_cache_source = folder;
                output.disk_cache_identity = stamp;
            end
        elseif strcmp(action,'save')
            if ~ffc_reconstruction_valid(scan) || ~isequaln(scan.disk_cache_identity,stamp), return; end
            if ~isfolder(cacheRoot), mkdir(cacheRoot); end
            cacheIdentity = stamp; %#ok<NASGU>
            cachedScan = scan; %#ok<NASGU>
            temp = [tempname(cacheRoot) '.mat'];
            cleanup = onCleanup(@() removeTemporary(temp)); %#ok<NASGU>
            save(temp,'cacheIdentity','cachedScan','-v7.3');
            if ~isequaln(stamp,identity(folder)), return; end
            movefile(temp,target,'f');
            trimCache(cacheRoot,target);
            output = true;
        else
            error('ffc_disk_cache:UnknownAction','Unknown cache action.');
        end
    catch ME
        warning('ffc_disk_cache:Unavailable','Reconstruction cache unavailable: %s',ME.message);
        output = [];
    end
end
function stamp = identity(folder)
    raw = dir(fullfile(folder,'data.dat'));
    header = fullfile(folder,'header.xml');
    stamp = [];
    if isempty(raw) || ~isfile(header), return; end
    serie = fullfile(folder,'Serie.xml');
    serieText = ''; if isfile(serie), serieText = fileread(serie); end
    stamp = struct('version',26,'folder',folder,'bytes',raw(1).bytes, ...
        'modified',raw(1).datenum,'header',fileread(header),'serie',serieText);
end
function removeTemporary(file)
    if isfile(file), delete(file); end
end
function trimCache(root,current)
% Bound cache growth, except a single file larger than the limit is retained.
    files = dir(fullfile(root,'*.mat'));
    files = files(~cellfun(@isempty,regexp({files.name},'^[0-9a-f]{64}\.mat$','once')));
    limit = getpref('FFCProcessing','DiskCacheLimitGiB',10)*1024^3;
    if ~isscalar(limit) || ~isfinite(limit) || limit <= 0, return; end
    total = sum([files.bytes]);
    [~,order] = sort([files.datenum]);
    for k = order
        if total <= limit, break; end
        file = fullfile(root,files(k).name);
        if strcmp(file,current), continue; end
        delete(file); total = total-files(k).bytes;
    end
end
