function [obj] = open_mrifile()
%Opens user-designated file, determines what system it came from and
%creates appropriate object containing the file data.

persistent lastFolder


try
if isempty(lastFolder)
    lastFolder = cd;
end
bakCD = cd;
cd(lastFolder);
catch
    cd('C:\') %default to C if all else fails
end


[nmrfilename,nmrpathname] = uigetfile({'*.MRD;*.mat;','MRI Files';...
    '*.*','All Files' },'Select MR Data File',...
    'MultiSelect', 'on');

if isequal(nmrfilename,0) || isequal(nmrpathname,0)
    return; % The user pressed cancel.
end

cd(bakCD);
if ischar(nmrpathname)
    lastFolder = nmrpathname;
end

if iscell(nmrfilename)
    for n=1:length(nmrfilename)
        [~,~,filetype] = bst_fileparts([nmrpathname nmrfilename{n}]);
        cd(nmrpathname)
        [fid,~] = fopen(fullfile(nmrfilename{n}));
        obj{n} = ImageReconCore(fid,filetype);
    end
else
    [~,~,filetype] = bst_fileparts([nmrpathname nmrfilename]);
    cd(nmrpathname)
    if strcmp(filetype,'.mat')
        file = load(nmrfilename);
        if strcmp(class(file.saveList{1}),'ImageReconCore')
            obj = file.saveList;
        else
            objects = cellfun(@class,file.saveList,'UniformOutput',false);
            whitelist = {'H9_se_nav_Elina','H9_se_multislice','h9_flash','h9_flash_v2','H9_se_nav_v6','H9_ge_lock','H9_se_propeller','H9_se_nav_v7','H9_se_nav_v8','H9_se_nav_v9','H9_se_nav_v10','H9_se_nav_fermi_v1','H9_se_nav_v7','H9_ir_se','H9_ir_se_nav_v4','H9_se_nav_v4'}; %these are imaging sequences, we discard all non-imaging sequences
            toprocess = find(ismember(objects,whitelist));
            indx = 0;
            for n=toprocess
                
                indx =indx+1;
                try
                    [fid,~] = fopen(fullfile(nmrfilename));
                    obj{indx} = ImageReconCore(fid,filetype,n);
                catch ME
                    disp(['Error while processing sequence number ' num2str(n) ':'])
                    disp(ME)
                end
            end
        end
    end
end

