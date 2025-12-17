classdef FIDclass
    %Imports sequence parameters and raw nmr data from SMIS .mrd file.
    %Performs necessary data reordering to account for multi-slice/multi-echo
    %sequences
    %Performs 2DFT on raw data.
    %
    %
    %Written by James Ross
    
    
    properties (SetAccess = public)
        sequence %the pulse sequence ppl name (eg gechoa, fse etc)
        nmrdata %raw k-space data
        AR
        AI
        phase
        image %unmodified image
        trueimage
        nfids %number of nmr excitations used to generate the file
        nmrdatatype %data type
        nmrechoes %number of echoes. Defaults to 1 for sequences if unused.
        nmrexperiments %number of experiments
        nmrsamples %number of samples in the frequency encoding direction
        nmrslices %number of slices
        nmrviews %number of phase encoding steps (2D)
        nmrviews2 %number of phase encoding steps (3D)
        times
        fields
        fid %file ID
        sortdata
        experiments
        echo
        PM
        param = struct %other sequence parameters go here, eg slice thickness etc.
        
    end
    
    methods
        function obj = FIDclass(fid) %Constructor function. Takes file ID (fid) and fills properties with the appropriate data from the file.
            % Read file into variables
            obj.fid = fid;
            % First, the header.
            obj.nmrsamples = double(fread(fid, 1, '*int32'));
            obj.nmrviews = double(fread(fid, 1, '*int32'));
            obj.nmrviews2 = fread(fid, 1, '*int32');
            obj.nmrslices = fread(fid, 1, '*int32');
            fseek(fid, 2, 'cof');
            obj.nmrdatatype = fread(fid, 1, '*int16'); % Yes, but assume complex floats.
            fseek(fid, hex2dec('98'), 'bof');
            obj.nmrechoes = fread(fid, 1, '*int32');
            obj.nmrexperiments = fread(fid, 1, '*int32');
            %if nviews ~= nsamples, error('This script requires a square imaging matrix'), end
            obj.nfids = double(obj.nmrviews*obj.nmrviews2*obj.nmrslices*obj.nmrechoes*obj.nmrexperiments);
            
            % Next, the data.
            fseek(fid, hex2dec('200'), 'bof');
            [AR, countr] = fread(fid, [obj.nmrsamples,obj.nfids], '*float32', 4);
            fseek(fid, -8*obj.nmrsamples*obj.nfids + 4, 'cof');
            [AI, counti] = fread(fid, [obj.nmrsamples,obj.nfids], '*float32', 4);
%             if countr ~= counti, error('Data length mismatch during read'), end
%             if countr ~= obj.nmrsamples*obj.nfids, error('Data length mismatch during read'), end
            obj.AR = double(AR);
            obj.AI = double(AI);
            obj.nmrdata = double(AR + AI*1i);
      
            % Finally, the parameters, and close the file.
            pb = fread(fid)';
            ps = char(pb);
            ps((pb==44)|(pb==10)|(pb==13)|(pb==0)) = 32;    % erase special characters
            ps = strtrim(regexp(ps,' :','split'));  % clean up whitespace
            P = regexp(ps,'\s+','split');
            for j=2:size(P,2)-1,
                if isequal(P{1,j}{1},'VAR') || isequal(P{1,j}{1},'VAR_ARRAY'), P{1,j}(1) = []; end %builds a params structure
                eval(['params.' genvarname(P{1,j}{1}) ' = P{1,j}(2:end);']);
            end
            obj.param =  params; %fills param property with the sequence parameters
            fclose(fid);
            [name,file] = fileparts(obj.param.PPL{:}); %decomposes pulse sequence path and file name into parts
            

            
            obj.sequence = file;
            
            
            
        end
    end
end

