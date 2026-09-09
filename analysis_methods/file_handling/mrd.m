function [MRDObject,name] = mrd(varargin)
%Opens user-designated file, determines what system it came from and
%creates appropriate object containing the file data.
numvarargs = length(varargin);
[nmrfilename,nmrpathname] = uigetfile2({'*.MRD;*.SDF','MRI Files';...
    '*.*','All Files' },'Select MR Data File',...
    'MultiSelect', 'on');
if iscell(nmrfilename)
    for i=1:length(nmrfilename)
        [fid, message] = fopen(fullfile(nmrpathname,nmrfilename{i}));
        
        if isequal(nmrfilename{i},0) || isequal(nmrpathname,0),
            return; % The user pressed cancel.
            
        end
        
        % Open file.
        % else
        %     optargs(1:numvarargs) = varargin;
        %             [nmrfilename] = optargs{:};%if the user doesn't enter a path, let them pick one
        %             [pathstr] = fileparts(nmrfilename);
        % cd(pathstr)
        % [fid, message] = fopen(fullfile(nmrfilename));
        % if fid == -1, error(message), end
        
        
        
        
        
        %parse filename into components
        [pathstr, name, ext] = fileparts(nmrfilename{i});
        %use file extension to determine what system the data is from
        if  strcmpi(ext,'.mrd'); %MRD files are from C10
            initMRDObject = C10Main(fid); %if MRD file create a C10 object
            %   newobject = feval(MRDObject.sequence, fid);
            [fid, message] = fopen(fullfile(nmrpathname,nmrfilename{i}));
            MRDObject{i} = feval(initMRDObject.sequence, fid); %call the class file appropriate to the sequence (eg for a Gechoa image, call class GECHOA)
        end
    end
    
    
    %         tempobj = MRDObject{1};
    %         for p=2:length(MRDObject)
    %             tempobj.rawdata = tempobj.rawdata+MRDObject{p}.rawdata;
    %         end
    %         tempobj.rawdata=tempobj.rawdata./length(MRDObject);
    %         MRDObject = tempobj;
else
    
    if numvarargs < 1
        [fid, message] = fopen(fullfile(nmrpathname,nmrfilename));
        
        if isequal(nmrfilename,0) || isequal(nmrpathname,0),
            return; % The user pressed cancel.
            
        end
        
        % Open file.
    else
        optargs(1:numvarargs) = varargin;
        [nmrfilename] = optargs{:};%if the user doesn't enter a path, let them pick one
        [pathstr] = fileparts(nmrfilename);
        cd(pathstr)
        [fid, message] = fopen(fullfile(nmrfilename));
        if fid == -1, error(message), end
    end
    
    %parse filename into components
    [pathstr, name, ext] = fileparts(nmrfilename);
    %use file extension to determine what system the data is from
    if  strcmpi(ext,'.mrd'); %MRD files are from C10
        initMRDObject = C10Main(fid); %if MRD file create a C10 object
        %   newobject = feval(MRDObject.sequence, fid);
        [fid, message] = fopen(fullfile(nmrpathname,nmrfilename));
        %         try
        MRDObject = feval(initMRDObject.sequence, fid); %call the class file appropriate to the sequence (eg for a Gechoa image, call class GECHOA)#
        %         catch exception
        %             MRDObject = feval('FCGE', fid);
        %     end
    end
    
    
end
