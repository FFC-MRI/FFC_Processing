function [obj] = open_datafile(fileName, filePath)
%Opens user-designated file, determines what system it came from and
%creates appropriate object containing the file data.

persistent previousFolder


try
if isempty(previousFolder)
    previousFolder = cd;
end
bakCD = cd;
cd(previousFolder);
catch
    cd('C:\') %default to C if all else fails
end


if isequal(fileName,0) || isequal(filePath,0)
    return; % The user pressed cancel.
end

cd(bakCD);
if ischar(filePath)
    lastFolder = filePath;
end

if iscell(fileName)
    for n=1:length(fileName)
        [~,~,filetype] = bst_fileparts([filePath fileName{n}]);
        cd(filePath)
        [fid,~] = fopen(fullfile(fileName{n}));
        obj{n} = ImageReconCore(fid,filetype);
    end
else
    [~,~,filetype] = bst_fileparts([filePath fileName]);
    cd(filePath)
    if strcmp(filetype,'.mat')
        file = load(fileName);
        if strcmp(class(file.saveList{1}),'ImageReconCore')
            obj = file.saveList;
        else
            objects = cellfun(@class,file.saveList,'UniformOutput',false);
            whitelist = {'H9_se_multislice_cardiac','H9_se_nav_v10_cardiac','H9_se_nav_v9_cardiac','H9_ir_multislice_se','H9_se_nav_Elina','H9_se_multislice','h9_flash','h9_flash_v2','H9_se_nav_v6','H9_ge_lock','H9_se_propeller','H9_se_nav_v7','H9_se_nav_v8','H9_se_nav_v9','H9_se_nav_v10','H9_se_nav_fermi_v1','H9_se_nav_v7','H9_ir_se','H9_ir_se_nav_v4','H9_se_nav_v4','H9_se_nav_longT1','H9_se_nav_v2','H9_ir_se_nav_v2'}; %these are imaging sequences, we discard all non-imaging sequences
            toprocess = find(ismember(objects,whitelist));
            indx = 0;
            for n=toprocess
                
                indx =indx+1;
                try
                    [fid,~] = fopen(fullfile(fileName));
                    obj{indx} = ImageReconCore(fid,filetype,n);
                catch ME
                    disp(['Error while processing sequence number ' num2str(n) ':'])
                    disp(ME)
                end
            end
        end
    end
end

