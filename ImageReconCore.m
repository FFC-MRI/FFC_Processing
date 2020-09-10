classdef ImageReconCore
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        sequence %the pulse sequence ppl name (eg gechoa, fse etc)
        displayname %user friendly sequence name eg Spin-Echo
        name %filename
        filelocation %fileID
        fileindex %which file is this in the savelist?
        rawdata
        originalcomplexkspace
        complexkspace %raw MRI data before reconstruction
        compleximage %complex images
        magimage
        magkspace
        scaledimages
        phaseimage
        nfids %number of nmr excitations used to generate the file
        nmrdatatype %data type
        echoes %number of echoes. Defaults to 1 for sequences if unused.
        experiments %number of experiments
        samples %number of samples in the frequency encoding direction
        slices %number of slices
        views %number of phase encoding steps (2D)
        views2 %number of phase encoding steps (3D)
        bladeangles %propellor blade angles
        twopointT1 %two pt acquisition strategy
        n_timepoints %number of evolution times
        n_fieldpoints %number of evolution fields
        timepoints %list of evolution times in seconds
        fieldpoints %list of evolution fields in mT
        thk %slice thickness
        fov %fov in mm
        TE %Echo time
        FlipAngle %indegrees
        averages %signal averages
        orientation %scan orientation
        phase_direction %0 vertical or 1 horizontal
        resolution_inplane %in plane resolution in mm
        resolution_throughplane %through plane resolution in mm
        bandwidth %acquisition bandwidth
        partialkspace %is this a half Fourier acquisition?
        T1Maps %Put processed T1/R1 maps here
        R1Maps %Put processed T1/R1 maps here
        window_function %windowing function
        window_size     %size of window
        fft_size        %size of 2dfft (symmetrical)
        denoise_filter  %denoising filter
        denoise_params  %filter kernel
        mask
        dispersioncurve
        R1T1 %R1 or T1
        fid %file ID
        backgroundselect %should the user define the background?
        reprocess %should the file be reprocessed?
        checkfit %should the user check the t1 fit?
        n_receivers %number of receivers
        multichannel_recon %how to deal with multiple channels?
        rot_theta %describes rotation applied using GUI
        param = struct %other sequence parameters go here       
    end
    
    
    methods
        function obj = ImageReconCore(fid,filetype,varargin) %Constructor function. Takes file ID (fid) and fills properties with the appropriate data from the file.
            % patch for full path argument (LB 25/06/12, 3 lines)
            switch filetype
                %% 
                case '.mat'
                    filename = fopen(fid);
                    obj.filelocation = filename;
                    file = load(filename);
                    if size(varargin)>0
                        index =varargin{1};
                    else
                        index = length(file.saveList);
                    end
                    file = file.saveList{index};
                    obj.fileindex = index;
                    obj.samples = cell2mat(file.pprParamList(strcmp('NO_SAMPLES', file.pprParamList(:,1)),3));
                    obj.views = cell2mat(file.pprParamList(strcmp('NO_VIEWS', file.pprParamList(:,1)),3));
                    obj.views2 = size(file.data,3);
                    obj.slices = cell2mat(file.pprParamList(strcmp('NO_SLICES', file.pprParamList(:,1)),3));
                    obj.echoes = size(file.data,5);
                    obj.experiments = size(file.data,6)*size(file.data,7);
                    obj.averages = size(file.data,8);
                    obj.n_receivers = size(file.data,9);
                  
                    
                    truedim1 = size(file.data,1); %try to deal with aborted or corrupted scans by looking at what we have vs what we expected
                    truedim2 = size(file.data,2);
                    if truedim1 ~=  obj.samples
                        file.data = padarray(file.data,abs(obj.samples-truedim1),0,'post');
                    end
                    if truedim2 ~=  obj.views
                        file.data = padarray(file.data,[0,abs(obj.views-truedim2)],0,'post');
                    end

                    obj.rawdata = reshape(file.data,obj.samples,obj.views,[]); %reshaping to the correct dimensions will happen later, for now stick with [2D,everything else] for preprocessing
                    obj.n_timepoints = size(file.data,6);
                    obj.n_fieldpoints = size(file.data,7);
                    obj.timepoints = file.waveformProfile.generator.Tevo.*1000;
                    obj.fieldpoints = file.waveformProfile.generator.Bevo.*1000;
                    obj.fov = cell2mat(file.pprParamList(strcmp('FOV', file.pprParamList(:,1)),2));
                    obj.TE = cell2mat(file.pprParamList(strcmp('te', file.pprParamList(:,2)),3));
                    obj.FlipAngle = cell2mat(file.pprParamList(strcmp('alpha', file.pprParamList(:,2)),3));
                    temp1 = cell2mat(file.pprParamList(strcmp('SLICE_THICKNESS', file.pprParamList(:,1)),3));
                    obj.thk = temp1(2);
                    temp1 = file.pprParamList(strcmp('SAMPLE_PERIOD', file.pprParamList(:,1)),3);
                    temp1 = temp1{1,1}{1,1};
                    obj.bandwidth = 10000/temp1;
                    phaseoverundersample = cell2mat(file.pprParamList(strcmp('oversample2', file.pprParamList(:,2)),3));
                    if obj.views ~= obj.samples && phaseoverundersample == 0
                        obj.partialkspace = 1;
                    else
                        obj.partialkspace = 0;
                    end
                    try
                        obj.bladeangles = cell2mat(file.pprParamList(strcmp('prop_angle', file.pprParamList(:,2)),3));
                    catch
                    end
                    temp1 = cell2mat(file.pprParamList(strcmp('X_ANGLE', file.pprParamList(:,1)),3));
                    xangle = temp1(2)./10;
                    temp1 = cell2mat(file.pprParamList(strcmp('Y_ANGLE', file.pprParamList(:,1)),3));
                    yangle = temp1(2)./10;
                    temp1 = cell2mat(file.pprParamList(strcmp('Z_ANGLE', file.pprParamList(:,1)),3));
                    zangle = temp1(2)./10;
                    obj.orientation = [xangle yangle zangle];
                    obj.phase_direction = cell2mat(file.pprParamList(strcmp('PHASE_ORIENTATION', file.pprParamList(:,1)),3));
                    obj.resolution_inplane = obj.fov./obj.samples;
                    obj.resolution_throughplane = obj.thk;
                    obj.twopointT1 = 0;
                    obj.window_function = 'none';
                    obj.window_size = 1;
                    obj.fft_size = obj.samples;
                    obj.denoise_filter = 'none';
                    obj.denoise_params = []; %filter kernel
                    obj.param = file.pprParamList;
                    [~,obj.sequence] = bst_fileparts(cell2mat(file.pprParamList(strcmp('PPL', file.pprParamList(:,1)),2)));
                    switch obj.sequence
                        case 'H9_ir_se_nav_v2'
                            obj.sequence = 'H9_ir_se';
                        case 'H9_ir_se_nav_v3'
                            obj.sequence = 'H9_ir_se';
                        case 'H9_ir_se_nav_v4'
                            obj.sequence = 'H9_ir_se';
                        case 'H9_se_multislice'
                        case 'H9_ge_looklocker_nav'
                            obj.sequence = 'H9_ge_looklocker';
                            
                    end
                    fclose(fid);
                    %% 
                case '.MRD'
                    %% 
                    if ischar(fid)
                        fid = fopen(fid);
                    end
                    % Read file into variables
                    obj.fid = fid;
                    % First, the header.
                    obj.samples = single(fread(fid, 1, '*int32'));
                    obj.views = single(fread(fid, 1, '*int32'));
                    obj.views2 = single(fread(fid, 1, '*int32'));
                    obj.slices = single(fread(fid, 1, '*int32'));
                    fseek(fid, 2, 'cof');
                    obj.nmrdatatype = fread(fid, 1, '*int32'); % Yes, but assume complex floats.
                    fseek(fid, hex2dec('98'), 'bof');
                    obj.echoes = single(fread(fid, 1, '*int32'));
                    obj.experiments = single(fread(fid, 1, '*int32'));%%-8;%%%%%%%%%%%%%%%%%%%%%%
                    obj.nfids = single(obj.views*obj.views2*obj.slices*obj.echoes*obj.experiments);
                    
                    % Next, the data.
                    fseek(fid, hex2dec('200'), 'bof');
                    [AR, ~] = fread(fid, [obj.samples,obj.nfids], '*float32', 4);
                    fseek(fid, -8*obj.samples*obj.nfids + 4, 'cof');
                    [AI, ~] = fread(fid, [obj.samples,obj.nfids], '*float32', 4);
                    %            if countr ~= counti, error('Data length mismatch during read'), end
                    %            if countr ~= obj.samples*obj.nfids, error('Data length mismatch during read'), end
                    AR = single(AR);
                    AI = single(AI);
                    obj.rawdata = (AR) + (AI)*1i;
                    
                    % Finally, the parameters, and close the file.
                    pb = fread(fid)';
                    ps = char(pb);
                    ps((pb==44)|(pb==10)|(pb==13)|(pb==0)) = 32;    % erase special characters
                    ps = strtrim(regexp(ps,' :','split'));  % clean up whitespace
                    P = regexp(ps,'\s+','split');
                    for j=2:size(P,2)-1
                        if isequal(P{1,j}{1},'VAR') || isequal(P{1,j}{1},'VAR_ARRAY'), P{1,j}(1) = []; end %builds a params structure
                        %                 eval(['params.' genvarname(P{1,j}{1}) ' = P{1,j}(2:end);']);
                        eval(['params.' matlab.lang.makeValidName(P{1,j}{1}) ' = P{1,j}(2:end);']);
                    end
                    obj.param =  params; %fills param property with the sequence parameters
                    
                    if isfield(obj.param,'two_pt_switch')
                        obj.twopointT1 = str2double(obj.param.two_pt_switch{1});
                    else
                        obj.twopointT1 =0;
                    end
                    
                    switch obj.twopointT1
                        case 1
                            if isfield(obj.param,'b_evol')
                                obj.fieldpoints = (nonzeros(cellfun(@str2num, obj.param.b_evol(2:end))));
                                obj.fieldpoints =  obj.fieldpoints(1:obj.experiments);
                            else
                                obj.fieldpoints = 200;
                            end
                            if isfield(obj.param,'t_evol')
                                obj.timepoints = (nonzeros(cellfun(@str2num, obj.param.t_evol(2:end))));
                                obj.timepoints =  obj.timepoints(1:obj.experiments)';
                                if isempty(obj.timepoints)
                                    obj.timepoints = 0;
                                else
                                    
                                end
                            else
                                obj.timepoints = 0;
                            end
                        otherwise
                            if isfield(obj.param,'b_evol')
                                obj.fieldpoints = unique(nonzeros(cellfun(@str2num, obj.param.b_evol(2:end))));
                                temp =length(obj.fieldpoints);
                                obj.fieldpoints = cellfun(@str2num, obj.param.b_evol(2:temp+1));
                            else
                                obj.fieldpoints = 200;
                            end
                            if isfield(obj.param,'t_evol')
                                obj.timepoints = flipud(unique(nonzeros(cellfun(@str2num, obj.param.t_evol(2:end)))));
                                if isempty(obj.timepoints)
                                    obj.timepoints = 0;
                                end
                            else
                                obj.timepoints = 0;
                            end
                    end
                    obj.n_fieldpoints = length(obj.fieldpoints);
                    obj.n_timepoints = obj.experiments./obj.n_fieldpoints;
                    %              obj.timepoints = [0.252840000000000,0.165969000000000,0.108945000000000,0.0715140000000000,0.0469430000000000,0.0308140000000000,0.0202270000000000,0.161758000000000,0.106181000000000,0.0696990000000000,0.0457520000000000,0.0300320000000000,0.0197140000000000,0.0129410000000000,0.138245000000000,0.0907470000000000,0.0595680000000000,0.0391020000000000,0.0256670000000000,0.0168480000000000,0.0110600000000000,0.127527000000000,0.0837120000000000,0.0549500000000000,0.0360700000000000,0.0236770000000000,0.0155420000000000,0.0102020000000000,0.0957850000000000,0.0628750000000000,0.0412730000000000,0.0270920000000000,0.0177840000000000,0.0116740000000000,0.00766300000000000,0.0657500000000000,0.0431600000000000,0.0283310000000000,0.0185970000000000,0.0122070000000000,0.00801300000000000,0.00526000000000000,0.0567000000000000,0.0372190000000000,0.0244310000000000,0.0160370000000000,0.0105270000000000,0.00691000000000000,0.00453600000000000,0.0541100000000000,0.0355190000000000,0.0233150000000000,0.0153050000000000,0.0100460000000000,0.00659500000000000,0.00432900000000000,0.166870000000000,0.109537000000000,0.0719020000000000,0.0471980000000000,0.0309820000000000,0.0203370000000000,0.0133500000000000,0.161053000000000,0.105718000000000,0.0693950000000000,0.0455530000000000,0.0299020000000000,0.0196280000000000,0.0128840000000000,0.156797000000000,0.102925000000000,0.0675620000000000,0.0443490000000000,0.0291120000000000,0.0191090000000000,0.0125440000000000,0.153473000000000,0.100742000000000,0.0661290000000000,0.0434090000000000,0.0284940000000000,0.0187040000000000,0.0122780000000000,0.150690000000000,0.0989160000000000,0.0649300000000000,0.0426220000000000,0.0279780000000000,0.0183650000000000,0.0120550000000000,0.148222000000000,0.0972960000000000,0.0638670000000000,0.0419240000000000,0.0275200000000000,0.0180640000000000,0.0118580000000000,0.146005000000000,0.0958410000000000,0.0629120000000000,0.0412960000000000,0.0271080000000000,0.0177940000000000,0.0116800000000000,0.144090000000000,0.0945840000000000,0.0620860000000000,0.0407550000000000,0.0267520000000000,0.0175610000000000,0.0115270000000000,0.142487000000000,0.0935320000000000,0.0613960000000000,0.0403020000000000,0.0264550000000000,0.0173650000000000,0.0113990000000000,0.141193000000000,0.0926820000000000,0.0608380000000000,0.0399350000000000,0.0262140000000000,0.0172080000000000,0.0112950000000000]'.*1000;
                    obj.thk = str2double(obj.param.SLICE_THICKNESS{3});
                    obj.fov = str2double(obj.param.FOV{1});
                    obj.resolution_inplane = obj.fov./obj.samples;
                    obj.resolution_throughplane = obj.thk;
                    obj.fieldpoints = flipud(obj.fieldpoints);
                    obj.timepoints = reshape(obj.timepoints,obj.n_timepoints,[])';
                    clear obj.timepoints;
                    %  obj.timepoints=[0.454940000000000,0.196028000000000,0.0844660000000000,0.0363950000000000;0.359515000000000,0.154910000000000,0.0667490000000000,0.0287610000000000;0.243638000000000,0.104980000000000,0.0452350000000000,0.0194910000000000;0.157710000000000,0.0679550000000000,0.0292810000000000,0.0126170000000000;0.102085000000000,0.0439870000000000,0.0189530000000000,0.00816700000000000;0.0660800000000000,0.0284730000000000,0.0122690000000000,0.00528600000000000].*1000;
                    fclose(fid);
                    [~,file] = bst_fileparts(obj.param.PPL{:}); %decomposes pulse sequence path and file name into parts
                    obj.sequence = file;
                    
                    %default recon params
                    obj.window_function = 'None';
                    obj.window_size = 1;
                    obj.fft_size = obj.samples;
                    obj.denoise_filter = 'None';
                    obj.denoise_params = []; %filter kernel
                    
                    switch obj.sequence
                        case 'H9_ir_se_nav_v2'
                            obj.sequence = 'H9_ir_se';
                        case 'H9_ir_se_nav_v3'
                            obj.sequence = 'H9_ir_se';
                        case 'H9_ir_se_nav_v4'
                            obj.sequence = 'H9_ir_se';
                        case 'H9_se_multislice'
                        case 'H9_ge_looklocker_nav'
                            obj.sequence = 'H9_ge_looklocker';
                    end
            end
            %% 
        end %function
        
        function obj = rename_sequence(obj)
            switch lower(strcat(char(regexp(obj.sequence,'[a-zA-Z]','match'))'))
                case 'hsenavv'
                    obj.displayname = 'Field-Cycling Spin-Echo';
                case 'hsemultislice'
                    obj.displayname = 'Multislice Spin-Echo';
                case 'hflash'
                    obj.displayname = 'Turbo Flash';
                otherwise
                    obj.displayname = obj.sequence; %for sequences we don't recognise
            end
        end
            
        function obj = preprocessing(obj) %here we do things like reorder kspace or do phase corrections. This should ideally only be done once per data set
           
            if isempty(obj.originalcomplexkspace)||(obj.reprocess==1)
             temp1 = load(obj.filelocation);
                 matfile = temp1.saveList{obj.fileindex};
%                 clear obj.originalcomplexkspace obj.complexkspace
                obj.views = size(obj.rawdata,2);
                A = double(obj.rawdata);
               
                
                
                A = kspace_reorder(A,obj);  %if necessary account for non-sequential ordering eg centre out
                
                      
                
                A = reshape(A,[obj.samples,obj.views,obj.echoes,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);
                
                if obj.echoes>1
%                   for r=1:obj.n_receivers
%                       for f=1:obj.n_fieldpoints
%                           for t=1:obj.n_timepoints
%                             for s=1:obj.slices
%                                 for v=1:obj.views
%                                     
%                     offset= FindOffsetFromFid(A(:,v,2,s,t,f,r),1/(1000*obj.bandwidth),0)/64;
%                              mult = exp(-1i*2*pi/128*(0:128-1)*offset)';
%                               
%                              A(:,v,1,s,t,f,r)=A(:,v,1,s,t,f,r).*mult;
%                                 end
%                             end
%                           end
%                       end
%                   end
                  
                    A(:,:,2,:,:,:) = []; %throw out the navigator data
                else
                    A = squeeze(A);
                end
%                   A = removespikes(A);
                A = reshape(A,obj.samples,obj.views,[]);
                if obj.samples ~= obj.views && obj.partialkspace ==1
                     A = padarray(double(A),[0,double(obj.samples-obj.views)],0,'post');
                     if mean(A(:,1,1))==0
                         A(:,1,:) = A(:,2,:);
                     end
                    %             AA = centre_kspace(AA);
                    obj.views = obj.samples;
                    %PF recon if we need to
                    try
                        for n=1:size(A,3)
                            [~, A(:,:,n)] = pocs(A(:,:,n),10,0);
                        end
                    catch
                    end
                end
                A = reshape(A,[obj.samples,obj.views,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);

               [correctedkspace] = correct_phase(A,obj.backgroundselect,obj.n_receivers);
               
                obj.complexkspace = reshape(correctedkspace,[obj.samples,obj.views,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);
                obj.originalcomplexkspace =obj.complexkspace;
%                 matfile.dataProcessed = obj.complexkspace;
%                 temp1.saveList{obj.fileindex} = matfile;
%                 saveList = temp1.saveList;
%                 [directory,~] = bst_fileparts(obj.filelocation);
%                 processed_dir = fullfile(directory,'ProcessedData.mat');
%                 save(processed_dir,'saveList');
            else
%                 obj.complexkspace = matfile.dataProcessed;
                obj.complexkspace = reshape( obj.complexkspace,[obj.samples,size(obj.complexkspace,2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]); %enforce dimensionality
                obj.originalcomplexkspace =obj.complexkspace;
            end
        end
            
            
        function obj = buildimages(obj)
            
            %phase encoding vertical 
            
            upscale_factor_read = double(obj.fft_size-obj.samples)/2;
            upscale_factor_phase = double((obj.fft_size*(obj.views/obj.samples))-obj.views)/2;
           
            obj.complexkspace = windowkspace(obj.originalcomplexkspace,obj.window_size,obj.window_function); %perform the kspace windowing first
            for n=1:obj.n_receivers
                for s=1:obj.slices
            noise(n,s,:) = obj.originalcomplexkspace(1,:,s,1,1,n); %used in multicoil recon
                end
            end

            obj = correct_orientation(obj);     
% obj.complexkspace = mean(obj.complexkspace,6);
            obj.compleximage=fft2c(padarray(obj.complexkspace ,[round(upscale_factor_phase) round(upscale_factor_read)],0));
%             tempdims = size(obj.compleximage);
%             temp = reshape(obj.compleximage,tempdims(1),tempdims(2),[]);
%  
%             mat = [0.87,0 0;0 1 0;0 0 1];
%             tform_t = affine2d(mat);
%             for n=16:20
%             temp2(:,:,n)= imwarp(temp(:,:,n),tform_t);
%             temp(:,:,n) = padarray(temp2(:,:,n),[0,(128-size(temp2,2))/2,0],0,'both');
%             end
%             
%             temp = reshape(temp,tempdims);
%              obj.compleximage = temp;
            obj.magimage = combine_channels(obj.compleximage,noise,obj.multichannel_recon);

%                obj.magimage = abs( obj.compleximage(:,:,:,:,1));
%             obj.magimage = mean(obj.magimage,6); %average multichannel data for now
            obj.magkspace = abs(fft2c(obj.magimage));
            obj.magimage = ffc_mri_filter(obj.magimage,obj.denoise_filter,obj.denoise_params);
            
            obj.phaseimage = angle(obj.compleximage);
            
        end %here the conversion from k-space to image space is performed, including FFT and multi-coil reconstruction
        
        
        function obj = MapRelaxation(obj)
            %              r1t1mapping;
            dim = size(obj.magimage,1);
            fields = obj.fieldpoints; %list of evolution times in seconds
            times = obj.timepoints;
            n_fields = obj.n_fieldpoints;
            clear obj.T1Maps
            clear obj.R1Maps
            obj.T1Maps = zeros(dim,dim,n_fields);
            obj.R1Maps = zeros(dim,dim,n_fields);
            %             %             mask = roipoly;
            %             switch obj.twopointT1
            %                 case 1
            %                     sz = size(obj.compleximage);
            %                     %                     temp = reshape(obj.compleximage,[dim*dim,[],size(fields)]);
            %                     for n=1:size(fields)
            %                         temp = obj.compleximage(:,:,:,:,:,n);
            %                         signal(:,n) = (abs(temp(:)));
            %                     end
            %                     for j = 1:size(signal,1)
            %                         t1maptemp(j,:) =  TwoPointMultifieldMethod(signal(j,:),fields,200,0.92,times,1,15,2);
            %                     end
            %
            %                     obj.T1maps = reshape(t1maptemp,sz);
            %
            %                 otherwise
            %                     clear obj.T1maps
            %                     image1 = abs(obj.compleximage(:,:,1,1,1));
            %                     dummymask = imbinarize(image1);
            %
            %                     [idx, centroids]=kmeans(double(image1(:)),5,'distance','sqEuclidean','Replicates',3);
            %                     segmented_images = cell(1,5);
            %                     for k = 1:5
            %                         color = zeros(size(image1));
            %                         color(idx==k) = image1(idx==k);
            %                         segmented_images{k} = color;
            %                     end
            %                     for n=1:5
            %                         t1mask(:,:,n) = imbinarize(segmented_images{n});
            %
            %                     end
            %                     for n=1:n_fields
            %                         tempt1mask = zeros(dim,dim);
            %                         for k=1:5
            %                             for t=1:size(times,2)
            %                                 temp = abs(obj.compleximage(:,:,:,t,n));
            %
            %                                 signaltemp(t,k,n) = mean(temp(t1mask(:,:,k)==1));
            %                             end
            %                         end
            mF = 0.05;
            maskFactor = mF;
            dims = size(obj.magimage);
            nbrow = size(obj.magimage,1);
            nbcol = size(obj.magimage,2);
            t1mask = zeros(nbrow, nbcol, 1);
            
            
            
            
            maskTmp = t1mask(:,:);
            maskTmp = medfilt2(maskTmp); % remove salt and pepper noise
            maskThreshold = maskFactor*max(max(abs(obj.magimage(:,:,1,1,1,1))));
            maskTmp(find(abs(obj.magimage(:,:,1,1,1,1))> maskThreshold)) = 1;
            t1mask = maskTmp;
            clear maskTmp
            imagestobeprocessed = obj.magimage;
            if obj.checkfit
                for n=1:n_fields
                t1map = multipointT1map(squeeze(imagestobeprocessed(:,:,:,:,n)),times(n,:),1,t1mask);
                T1Maps(:,:,1,n)=t1map(:,:,:,1);
                end
            else
            parfor n=1:n_fields
                t1map = multipointT1map(squeeze(imagestobeprocessed(:,:,:,:,n)),times(n,:),0,t1mask);
                T1Maps(:,:,1,n)=t1map(:,:,:,1);
            end
            end
            
            
            obj.T1Maps = T1Maps;
            obj.R1Maps = 1000./T1Maps;
            
            %                         %                    x = times(1,:);
            %                         for k=1:5
            %                             y = signaltemp(:,k,n);
            %                             x = times(n,:);
            %                             [~,index] = min(y);
            %                             if index~=1
            %                                 y(1:index-1)=-1.*y(1:index-1);
            %                                 y = [y; -y(index)];
            %                                 x = [x x(index)];
            %                                 %                                 obj.T1maps(:,:,1,n) = T1map(squeeze(obj.compleximage(:,:,:,:,n)),times(1,:));
            %                             end
            %
            %                             signal = y;
            %                             invtimes = x./1000;
            %
            %                             [xData, yData] = prepareCurveData(double(invtimes(:)), double(signal(:)));
            %
            %                             % Set up fittype and options.
            %                             % if n==1
            %                             ft = fittype( '((M0-Minf)*exp(-x*R1)+Minf)', 'independent', 'x', 'dependent', 'y' );
            %                             opts = fitoptions( 'Method', 'NonlinearLeastSquares','Robust','on');
            %                             opts.Display = 'Off';
            %                             opts.Lower = [-Inf,-Inf,0];
            %                             opts.Upper = [Inf,Inf,Inf];
            %
            %                             [xs,ord] = sort(xData);
            %                             ys = yData(ord);  % obtain the signed magnitude
            %
            %                             T1est = 0.1;
            %                             opts.Startpoint = [ys(1),ys(end),1/T1est];
            %                             [fitob,gof] = fit(xs(:),ys(:),ft,opts);
            %
            %
            %                             gofmaptemp(n) = gof.rsquare;
            %                             error = confint(fitob,0.63);
            %                             errors(n) = diff(error(:,2));
            %                             R1dispersion(k,n)= fitob.R1;
            %
            %                         end
            %
            %                     end
            %                     for n=1:n_fields
            %                         for k=1:5
            %                             tempt1mask(t1mask(:,:,k)==1) =  1000./R1dispersion(k,n);%t1 instead of r1
            %                         end
            %                         obj.T1maps(:,:,n)=tempt1mask;
            %
            %                     end
            %                     obj.T1maps = obj.T1maps.*dummymask;
            %
            %             end
            
        end %image based T1 analysis methods are contained here
        
        function obj = T1dispersion(obj)
            fields = obj.fieldpoints; %list of evolution times in seconds
            times = obj.timepoints;
            n_fields = obj.n_fieldpoints;
            R1dispersion = zeros(1,n_fields);
              signal = abs(((squeeze(obj.compleximage(:,:,:,:,:)))));
            signal = signal.*repmat(obj.mask,[1 1 size(signal,3)]);
                 signal = sum(sum(signal,1),2);
                 signal = squeeze(signal);
%              [fitobject,gof,R1out] = fit_relaxation(signal,times./1000,fields,0.200);
          
            for n=1:n_fields
                invtimes = times(n,:);
                signal = abs(((squeeze(obj.compleximage(:,:,:,:,n)))));
                
                signal = signal.*repmat(obj.mask,[1 1 size(signal,3)]);
                signal = sum(sum(signal,1),2);
                signal = squeeze(signal);
                x = invtimes;
                y = signal;
                
                [~,index] = min(y);
                if index~=1
                    y(1:index-1)=-1.*y(1:index-1);
                    y = [y; -y(index)];
                    x = [x x(index)];
                end
                %
                signal = y(1:end);
                invtimes = x(1:end)./1000;
                
                [xData, yData] = prepareCurveData(double(invtimes(:)), double(signal(:)));
                
                % Set up fittype and options.
                % if n==1
                ft = fittype( '((M0-Minf)*exp(-x*R1)+Minf)', 'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','Robust','on');
                opts.Display = 'Off';
                opts.Lower = [-Inf,-Inf,0];
                opts.Upper = [Inf,Inf,Inf];
                
                % else
                % ft = fittype( 'a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
                % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                % opts.Display = 'Off';
                % opts.Lower = [-Inf 10];
                % opts.StartPoint = [-10000 10 ];
                % opts.Upper = [Inf 1000];
                % end
                %
                
                [xs,ord] = sort(xData);
                ys = yData(ord);  % obtain the signed magnitude
                if n>1
                    T1est = 1./R1dispersion(n-1);
                else
                    T1est = 0.2;
                end
                %                 Weights = yData;
                %                 opts.Weights = Weights;
                opts.Startpoint = [ys(1),ys(end),1/T1est];
                [fitob,gof] = fit(xs(:),ys(:),ft,opts);
                
                
                gofmaptemp(n) = gof.rsquare;
                error = confint(fitob,0.63);
                errors(n) = diff(error(:,2));
                
                
                R1dispersion(n)= fitob.R1;
                
                gofmaptemp(n) = gof.rsquare;
                error = confint(fitob,0.63);
                errors(n) = diff(error(:,2));
                if obj.checkfit ==1
                    h=figure;
                    
                    plot(fitob,xs(:),ys(:));
                    ylabel('Signal (AU)')
                    xlabel('Evolution Time (s)');
                    try
                        w = waitforbuttonpress;
                        close(h);
                    catch
                    end
                    
                end
                
            end
            T1dispersion = 1./R1dispersion.*1000; %in ms
            fieldHz = fields.*gammaH*1e-9;
            fields = fields./1e3;
%               figure, scatter(fieldHz, R1out);
            if obj.R1T1 ==1
                disp(['Fields (MHz): ' num2str(fieldHz')])
                disp(['R1 (s^-1): ' num2str(R1dispersion)])
                figure,scatter(fieldHz,R1dispersion,'o');
                set(gca,'xscale','log','yscale','log');
                xlabel('Evolution Field (MHz)')
                ylabel('R_1 (s^-^1)');
            else
                disp(['Fields (T): ' num2str(fields')])
                disp(['T1 (ms): ' num2str(T1dispersion)])
                figure,scatter(fields,T1dispersion,'o');
                set(gca,'xscale','log','yscale','log');
                xlabel('Evolution Field (T)')
                ylabel('T_1 (ms)');
                obj.dispersioncurve = R1dispersion;
            end
        end %ROI based T1 methods are contained here
    end
    
    
end


