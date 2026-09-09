[nmrfilename,nmrpathname] = uigetfile({'*.MRD;','MRI Files';...
    '*.*','All Files' },'Select MR Data File','MultiSelect', 'on');
if isequal(nmrfilename,0) || isequal(nmrpathname,0),
    return; % The user pressed cancel.
end
cd(nmrpathname)
if iscell(nmrfilename) %do we want to open many files?
    
    
    for i = 1:length(nmrfilename)
        
        [fid, ~] = fopen(fullfile(nmrpathname,nmrfilename{i}));
        if fid == -1, error(message), end
        
        [fid, message] = fopen(fullfile(nmrpathname,nmrfilename{i}));%open file
        [pathstr, name, ext] = fileparts(nmrfilename{i});
        v = genvarname([num2str(i) 'G FID']);   %generate a unique name for the file
        object = FIDclass(fid);                 %build an object to put the data in
        FID.(genvarname(name)) = object;    
      
        
        
    end
    
else
    [fid, message] = fopen(fullfile(nmrpathname,nmrfilename));
    [pathstr, name, ext] = fileparts(nmrfilename);
    v = genvarname([name 'G phase']);
    object = FIDclass(fid);
    %sort data by experiment/echoes JR 20/05/13
    A = object.nmrdata;         %raw complex data
    S = size(object.nmrdata,2)/object.nmrexperiments;   
    E = object.nmrechoes;   %number of echoes
    R = S/object.nmrechoes; %lines per echo
    try
        object.sortdata=reshape(A,[size(A,1),object.nmrechoes,object.nmrexperiments]);
    catch
        try
            object.sortdata=reshape(A,[size(A,1),object.nmrechoes,object.nmrviews,object.nmrexperiments]);
        catch
            object.sortdata=reshape(A,[size(A,1),object.nmrviews,object.nmrviews2,object.nmrslices,object.nmrexperiments]);
        end
    end
    
    FID.(genvarname(name)) = object;


end

clear name message file fid ext nmrpathname nmrfilename v pathstr


