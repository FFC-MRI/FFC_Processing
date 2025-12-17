function obj = open_mrifile()
% Open MRI data for reconstruction.
% - Option 1: choose a folder -> reconstruct all .dat files (recursive), excluding any under a 'proc' folder
% - Option 2: cancel folder dialog -> choose file(s) via uigetfile (supports .MRD, .mat, .dat)
%
% Returns:
%   obj : cell array of ImageReconCore objects (or empty {} on cancel / no valid files)

persistent lastFolder
obj = {};

% --- Initialise lastFolder safely ---
try
    if isempty(lastFolder) || ~isfolder(lastFolder)
        lastFolder = pwd;
    end
catch
    lastFolder = pwd;
end

%% =========================
%  1) Folder mode (batch .dat)
%  =========================
rootFolder = uigetdir(lastFolder, 'Select folder containing .dat files (Cancel to select files)');
if ~isequal(rootFolder, 0)
    if isstring(rootFolder); rootFolder = char(rootFolder); end
    lastFolder = rootFolder;

    % Find .dat files recursively
    datList = dir(fullfile(rootFolder, '**', '*.dat'));

    % Exclude anything under a folder named exactly 'proc'
    if ~isempty(datList)
        isUnderProc = false(numel(datList), 1);
        for k = 1:numel(datList)
            parts = regexp(datList(k).folder, filesep, 'split');
            isUnderProc(k) = any(strcmpi(parts, 'proc'));
        end
        datList = datList(~isUnderProc);
    end

    if isempty(datList)
        warning('No .dat files found under: %s (excluding ''proc'' folders)', rootFolder);
        obj = {};
        return;
    end

    n = numel(datList);
    obj  = cell(1, n);
    ok   = false(1, n);
    errs = cell(1, n);

    % Start pool if needed (optional; comment out if you prefer manual control)
    if isempty(gcp('nocreate'))
        try
            parpool;
        catch
            % If Parallel Toolbox not available, fall back to serial loop
            warning('Parallel pool could not be started. Falling back to serial processing.');
            for k = 1:n
                try
                    thisDatFull   = fullfile(datList(k).folder, datList(k).name);
                    thisDatFolder = datList(k).folder;

                    [data, par, data_name] = read_raw_data_Cameleon(thisDatFolder, []);
                    file = struct();
                    file.data = data;
                    file.par  = par;
                    file.data_name = data_name; %#ok<NASGU>

                    obj{k} = ImageReconCore(file, '.dat', 1);
                    ok(k)  = true;
                catch ME
                    ok(k) = false;
                    errs{k} = sprintf('Failed to reconstruct %s\n  %s', thisDatFull, ME.message);
                    obj{k} = [];
                end
            end

            obj = obj(ok);
            bad = find(~ok);
            for i = 1:numel(bad)
                if ~isempty(errs{bad(i)})
                    warning('%s', errs{bad(i)});
                end
            end
            return;
        end
    end

    % PARALLEL loop
    parfor k = 1:n
        thisDatFull   = fullfile(datList(k).folder, datList(k).name);
        thisDatFolder = datList(k).folder;

        try
            % Important: read_raw_data_Cameleon expects a dataset FOLDER in your current usage
            [data, par, data_name] = read_raw_data_Cameleon(thisDatFolder, []);

            % Keep all variables local per worker
            file = struct();
            file.data = data;
            file.par  = par;
            file.data_name = data_name; %#ok<NASGU>

            obj{k} = ImageReconCore(file, '.dat', 1);
            ok(k)  = true;

        catch ME
            ok(k) = false;
            errs{k} = sprintf('Failed to reconstruct %s\n  %s', thisDatFull, ME.message);
            obj{k} = [];
        end
    end

    % Compact results and warn (on client) for failures
    obj = obj(ok);

    bad = find(~ok);
    for i = 1:numel(bad)
        if ~isempty(errs{bad(i)})
            warning('%s', errs{bad(i)});
        end
    end

    return;
end

%% =========================
%  2) File mode (original uigetfile)
%  =========================
bakCD = pwd;
try
    cd(lastFolder);
catch
    cd('C:\');
end

[nmrfilename, nmrpathname] = uigetfile({ ...
    '*.MRD;*.mat;*.dat','MRI Files'; ...
    '*.*','All Files' }, ...
    'Select MR Data File', ...
    'MultiSelect','on');

cd(bakCD);

if isequal(nmrfilename,0) || isequal(nmrpathname,0)
    obj = {};
    return;
end

if ischar(nmrpathname)
    lastFolder = nmrpathname;
end

% Normalise selection into a cell array of filenames
if ischar(nmrfilename)
    fileList = {nmrfilename};
else
    fileList = nmrfilename;
end

out = {};  %#ok<NASGU>
obj = {};
outIdx = 0;

for n = 1:numel(fileList)
    thisFile = fileList{n};
    fullpath = fullfile(nmrpathname, thisFile);

    [~,~,filetype] = bst_fileparts(fullpath);

    switch lower(filetype)

        case '.mat'
            file = load(fullpath);

            % If saveList contains ImageReconCore already, return it as a cell array
            if isfield(file, 'saveList') && ~isempty(file.saveList)
                if isa(file.saveList, 'ImageReconCore')
                    outIdx = outIdx + 1;
                    obj{outIdx} = file.saveList; %#ok<AGROW>
                elseif iscell(file.saveList)
                    % Filter to whitelist imaging sequences (as in your original code)
                    objects = cellfun(@class, file.saveList, 'UniformOutput', false);
                    whitelist = {'H9_ge_nav_v1','H9_se_multislice_cardiac','H9_se_nav_v10_cardiac', ...
                        'H9_se_nav_v9_cardiac','H9_se_nav_v9_t1rho','H9_ir_multislice_se','H9_se_nav_Elina', ...
                        'H9_se_multislice','h9_flash','h9_flash_v2','H9_se_nav_v6','H9_ge_lock','H9_se_propeller', ...
                        'H9_se_nav_v7','H9_se_nav_v8','H9_se_nav_v9','H9_se_nav_v10','H9_se_nav_fermi_v1', ...
                        'H9_se_nav_v7','H9_ir_se','H9_ir_se_nav_v4','H9_se_nav_v4','H9_se_nav_longT1', ...
                        'H9_se_nav_v2','H9_ir_se_nav_v2'};

                    toprocess = find(ismember(objects, whitelist));

                    for seqIdx = toprocess
                        try
                            fid = fopen(fullpath);
                            outIdx = outIdx + 1;
                            obj{outIdx} = ImageReconCore(fid, filetype, seqIdx); %#ok<AGROW>
                        catch ME
                            warning('Error while processing sequence number %d in %s:\n  %s', seqIdx, fullpath, ME.message);
                        end
                    end
                else
                    warning('Loaded .mat does not contain a recognised saveList format: %s', fullpath);
                end
            else
                warning('No saveList found in .mat file: %s', fullpath);
            end

        case '.dat'
            % For file picker, treat the containing folder as the dataset folder
            thisFolder = fileparts(fullpath);

            try
                [data, par, data_name] = read_raw_data_Cameleon(thisFolder, []);
                file = struct();
                file.data = data;
                file.par  = par;
                file.data_name = data_name; %#ok<NASGU>

                outIdx = outIdx + 1;
                obj{outIdx} = ImageReconCore(file, '.dat', 1); %#ok<AGROW>

            catch ME
                warning('Failed to reconstruct %s\n  %s', fullpath, ME.message);
            end

        otherwise
            % .MRD etc
            try
                fid = fopen(fullpath);
                outIdx = outIdx + 1;
                obj{outIdx} = ImageReconCore(fid, filetype); %#ok<AGROW>
            catch ME
                warning('Failed to open %s\n  %s', fullpath, ME.message);
            end
    end
end

end
