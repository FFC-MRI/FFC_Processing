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
        complexcombined
        magimage
        magkspace
        scaledimages
        phaseimage
        psirrealimage    % <-- use this for fitting (signed)
        psirimagimage                 % residuals (diagnostic)
        psirphi_ref                  % saved for QC
        psirweights
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
        multichannel_recon %what channels are selected
        rot_theta %describes rotation applied using GUI
        flipH %mirror horizontally
        flipV %mirror vertically
        param = struct %other sequence parameters go here
        SENSE %undersampled
        registration_tform %transformation used for image registration
        TwoDimensional
        recon2d
        geomT
        selected_coils_original
        biascorrect
        combination
        probe %the RF coil used
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
                    if isempty(cell2mat(file.pprParamList(strcmp('NO_VIEWS_2', file.pprParamList(:,1)),3)))
                        obj.views2 = 1;
                    else
                        obj.views2 = cell2mat(file.pprParamList(strcmp('NO_VIEWS_2', file.pprParamList(:,1)),3));
                    end

                    obj.slices = cell2mat(file.pprParamList(strcmp('NO_SLICES', file.pprParamList(:,1)),3));

                    if isempty(cell2mat(file.pprParamList(strcmp('NO_ECHOES', file.pprParamList(:,1)),3)))
                        obj.echoes = 1;
                    else
                        obj.echoes = cell2mat(file.pprParamList(strcmp('NO_ECHOES', file.pprParamList(:,1)),3));
                    end



                    obj.experiments = size(file.data,6)*size(file.data,7);
                    obj.averages = size(file.data,8);
                    obj.n_receivers = size(file.data,9);
                    obj.recon2d = obj.samples;

                    truedim1 = size(file.data,1); %try to deal with aborted or corrupted scans by looking at what we have vs what we expected
                    truedim2 = size(file.data,2);

                    %                     if truedim1 > obj.samples
                    %                         file.data = file.data(1:obj.samples,:,:,:,:,:,:,:); %deal with that strange 4x 0 padding issue the new evo produces.
                    %                     end
                    %
                    if truedim1 <  obj.samples
                        file.data = padarray(file.data,abs(obj.samples-truedim1),0,'post');
                    end
                    if truedim2 <  obj.views
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
                    obj.SENSE = 0;
                    try
                        obj.SENSE = cell2mat(file.pprParamList(strcmp('sense_onoff', file.pprParamList(:,2)),3));
                    catch
                    end
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
                    xangle = temp1(2:obj.slices+1)./10;
                    temp1 = cell2mat(file.pprParamList(strcmp('Y_ANGLE', file.pprParamList(:,1)),3));
                    yangle = temp1(2:obj.slices+1)./10;
                    temp1 = cell2mat(file.pprParamList(strcmp('Z_ANGLE', file.pprParamList(:,1)),3));
                    zangle = temp1(2:obj.slices+1)./10;
                    obj.orientation = [xangle;yangle;zangle];
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
                    obj.TwoDimensional =1;
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

                case '.dat'
                    obj.rawdata = fid.data;
                    obj.orientation = fid.par.cameleon.IMAGE_ORIENTATION_SUBJECT;
                    obj.samples = fid.par.md1d;
                    obj.views = fid.par.md2d;
                    obj.slices = fid.par.md3d;
                    obj.experiments =fid.par.md4d;
                    obj.n_receivers = fid.par.ncoils;
                    obj.echoes = 1;
                    obj.TE = fid.par.te*1000;
                    obj.FlipAngle = fid.par.fa;
                    obj.fov = fid.par.cameleon.FIELD_OF_VIEW*1000;
                    obj.TwoDimensional = fid.par.cameleon.MULTI_PLANAR_EXCITATION;
                    obj.thk = fid.par.cameleon.SLICE_THICKNESS*1000;
                    obj.resolution_inplane = fid.par.cameleon.RESOLUTION_FREQUENCY;
                    obj.resolution_throughplane = obj.thk;
                    obj.sequence = fid.par.serieName;
                    obj.partialkspace =0;
                    if fid.par.cameleon.USER_ZERO_FILLING_2D ~=0
                        obj.partialkspace =1;
                    end
                    obj.geomT = eye(3);

                    obj.recon2d = fid.par.cameleon.USER_MATRIX_DIMENSION_2D;

                    obj.averages =1;
                    obj.twopointT1 = 0;
                    obj.window_function = 'none';
                    obj.window_size = 1;
                    obj.fft_size = obj.samples;
                    obj.denoise_filter = 'none';
                    obj.denoise_params = []; %filter kernel
                    obj.param = fid.par.cameleon;

                    if fid.par.cameleon.MULTI_PLANAR_EXCITATION ==1
                        obj.TwoDimensional = 1;
                    else
                        obj.TwoDimensional = 0; %ie 3d
                    end



                    switch obj.param.TRANSFORM_PLUGIN

                        case 'Sequential4DCine'
                            obj.n_timepoints =  obj.param.CARDIAC_NB_OF_PHASE;
                            obj.timepoints = obj.param.ACQUISITION_TIME_OFFSET;
                            obj.n_fieldpoints =1;
                            obj.fieldpoints = 200;

                        otherwise
                            obj.n_timepoints=1;
                            obj.n_fieldpoints =1;
                            obj.timepoints = 1;
                            obj.fieldpoints = 200;    
                    end

                    if isfield(obj.param,'SE_TYPE')
                        if strcmp(obj.param.SE_TYPE,'MultiEcho')
                            obj.n_timepoints = obj.param.ECHO_TRAIN_LENGTH;
                            obj.timepoints = obj.param.ECHO_TIME;
                        end
                    end

                    if isfield(obj.param,'FCI_LIST_OF_BEVO')
                        obj.n_fieldpoints=size(obj.param.FCI_LIST_OF_BEVO,2);
                        obj.n_timepoints = obj.param.FCI_NUMBER_OF_TEVO;
                        obj.timepoints = [];
                        obj.fieldpoints = []; 
                        if obj.param.FCI_T1_LOG_DIST ==1
                            for n=1:size(obj.param.FCI_LIST_OF_BEVO,2)
                                obj.timepoints(n,:) = logspace(log10(obj.param.FCI_LIST_OF_T1_ESTIMATES(n).*obj.param.FCI_T1_FACTOR_LOWER),log10(obj.param.FCI_LIST_OF_T1_ESTIMATES(n).*obj.param.FCI_T1_FACTOR_UPPER),obj.param.FCI_NUMBER_OF_TEVO);
                            end
                        else
                            for n=1:size(obj.param.FCI_LIST_OF_BEVO,2)
                                obj.timepoints(n,:) = linspace(obj.param.FCI_LIST_OF_T1_ESTIMATES(n).*obj.param.FCI_T1_FACTOR_LOWER,obj.param.FCI_LIST_OF_T1_ESTIMATES(n).*obj.param.FCI_T1_FACTOR_UPPER,obj.param.FCI_NUMBER_OF_TEVO);
                            end
                        end
                        obj.timepoints = obj.timepoints.*1000;
                        obj.fieldpoints=obj.param.FCI_LIST_OF_BEVO;
                      
                    end
                obj.probe = obj.param.PROBES;
      
     
                  


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
                try
                    temp1 = load(obj.filelocation);
                    matfile = temp1.saveList{obj.fileindex};
                catch
                end
                obj.views = size(obj.rawdata,2);
                A = single(obj.rawdata);
                A = reshape(A,[obj.samples,obj.views,obj.echoes,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.averages,obj.n_receivers]);
                A = mean(A,7);
                if 1 ==1
                    A = kspace_reorder(A,obj);  %if necessary account for non-sequential ordering eg centre out

                    if obj.echoes>1
                        A(:,:,2,:,:,:) = []; %throw out the navigator data
                    else
                        A = squeeze(A);

                    end

                    if obj.recon2d == 0
                    end

                    %   A = removespikes(A);
                    A = reshape(A,obj.samples,obj.views,[]);

                    % if obj.samples ~= obj.views && obj.partialkspace ==1
                    %     try
                    %         A = padarray(single(A),[0,single(obj.recon2d-obj.views)],0,'post');
                    %         obj.views = obj.recon2d;
                    %     catch
                    %         A = padarray(single(A),[0,single(obj.samples-obj.views)],0,'post');
                    %         obj.views = obj.samples;
                    %     end
                    %
                    %     if mean(A(:,1,1))==0
                    %         A(:,1,:) = A(:,2,:); %for some reason the console occasionally fills the first line of kspace with 0s which messes pocs up. This fixes it.
                    %     end
                    %     %             AA = centre_kspace(AA);
                    %
                    %     %PF recon if we need to
                    %     % try
                    %     %     for n=1:size(A,3)
                    %     %         [~, A(:,:,n)] = pocs(A(:,:,n),10,0);
                    %     %     end
                    %     % catch
                    %     end

                    % --- Symmetric undersampling about k0 (NOT partial Fourier) ---
                    % Phase-encode dimension is dim 2

                    fullViews = double(obj.param.USER_MATRIX_DIMENSION_2D);        % nominal
                    acqViews  = double(obj.param.ACQUISITION_MATRIX_DIMENSION_2D); % collected

                    % If you trust obj.views more than the header, you can cross-check:
                    % acqViews = double(obj.views);

                    % Optional sanity check against USER_PARTIAL_PHASE (doesn't drive padding)
                    if isfield(obj.param,'USER_PARTIAL_PHASE')
                        pct = double(obj.param.USER_PARTIAL_PHASE);
                        % expectedAcq = round(fullViews * pct/100);
                        % (don't enforce; some consoles round differently)
                    end

                    curViews = size(A,2);

                    % If the array doesn't match the header, prefer the actual array size.
                    % (This avoids padding the wrong amount if obj.views/header got stale.)
                    if curViews ~= acqViews
                        acqViews = curViews;
                    end

                    targetViews = fullViews;

                    if acqViews < targetViews
                        miss = targetViews - acqViews;

                        pre  = floor(miss/2);
                        post = miss - pre;   % handles odd miss cleanly

                        A = single(A); % ensure consistent type
                        if pre > 0
                            A = padarray(A, [0 pre],  0, 'pre');
                        end
                        if post > 0
                            A = padarray(A, [0 post], 0, 'post');
                        end

                        obj.views = targetViews;

                    elseif acqViews > targetViews
                        % Rare case: acquired more than nominal; crop symmetrically about k0
                        extra = acqViews - targetViews;
                        pre  = floor(extra/2);
                        post = extra - pre;

                        A = A(:, (1+pre):(end-post), :, :, :, :);  % keeps all higher dims
                        obj.views = targetViews;

                    else
                        % already the right size
                        obj.views = acqViews;
                    end

           

                    A = reshape(A,[obj.samples,obj.views,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);
                    %  [correctedkspace] = correct_phase(A,obj.backgroundselect,obj.n_receivers);
                    correctedkspace = A;
                    % fprintf("kspace delta: %.3g\n", norm(correctedkspace(:)-A(:)) / norm(A(:)));
                else
                    correctedkspace = reshape(A,[obj.samples,obj.views,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);
                end
                obj.complexkspace = reshape(correctedkspace,[obj.samples,obj.views,obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);
                %so we can undo windowing etc keep an untouched version of kspace prior to FFT
            if strcmp(string(obj.probe), "Body Coil 4Rx Coil")
            obj.originalcomplexkspace = noise_whiten(obj.complexkspace);
                else
             obj.originalcomplexkspace = obj.complexkspace;
            end
               
            else
                obj.complexkspace = reshape( obj.complexkspace,[obj.samples,size(obj.complexkspace,2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]); %enforce dimensionality
                if strcmp(string(obj.probe), "Body Coil 4Rx Coil")
            obj.originalcomplexkspace = noise_whiten(obj.complexkspace);
                else
             obj.originalcomplexkspace = obj.complexkspace;
            end
            end

        end


        function obj = buildimages(obj)

            %BUILDIMAGES  Convert k-space to image space, perform coil combination, optional bias correction, filtering.
            %
            
            obj.compleximage =[];
            obj.complexkspace = obj.originalcomplexkspace;

            obj = correct_orientation(obj);
   
            % ----------- Precompute padding -----------
            upscale_factor_read  = double(obj.fft_size - obj.samples)/2;
            upscale_factor_phase = double((obj.fft_size*(obj.views/obj.samples)) - obj.views)/2;

            padPhase = max(0, round(upscale_factor_phase));
            padRead  = max(0, round(upscale_factor_read));
            
            % ----------- Window k-space first -----------
            obj.complexkspace = windowkspace(obj.complexkspace, obj.window_size, obj.window_function);

            % ----------- Pad k-space (phase/read are dims 1/2 in your convention) -----------
            kpad_full = padarray(obj.complexkspace, [padPhase padRead], 0);

            % ----------- Coil selection (trim once and renumber) -----------
            % obj.multichannel_recon may be:
            %   - scalar 1 (meaning "all coils")
            %   - logical mask over the ORIGINAL coil dimension
            %   - numeric index list into the ORIGINAL coil dimension (e.g. [2 3 4 5 6 7])
            opts_orig   = obj.multichannel_recon;      % numeric index list (e.g. [2 3 4 5 6 7])
            nCoils_orig = size(kpad_full, 6);
            probeIsQBC = false;
            try
                probeIsQBC = isfield(obj.probe) && strcmp(string(obj.probe), "Quadrature BirdCage");
            catch
                probeIsQBC = false;
            end

            opts_orig = obj.multichannel_recon;

            userDidNotSelect = 0;

            if probeIsQBC && userDidNotSelect
                if nCoils_orig >= 8
                    obj.multichannel_recon = 8;   % use only channel 8
                else
                    obj.multichannel_recon = 1;   % fallback to channel 1
                end
                opts_orig = obj.multichannel_recon; % update local copy
            end
            if isempty(opts_orig)
                sel_orig = 1:nCoils_orig;
            else
                if ~isnumeric(opts_orig)
                    error('buildimages:coilSelNotNumeric', ...
                        'multichannel_recon must be a numeric index list. Got %s.', class(opts_orig));
                end

                sel_orig = unique(round(opts_orig(:).'));   % row vector, unique
                sel_orig(isnan(sel_orig)) = [];

                if isempty(sel_orig)
                    error('buildimages:coilSelEmptyAfterParse', ...
                        'multichannel_recon became empty after parsing.');
                end

                if any(sel_orig < 1) || any(sel_orig > nCoils_orig)
                    error('buildimages:coilIndexOutOfRange', ...
                        'Selected coil indices [%s] exceed available coils 1..%d.', num2str(sel_orig), nCoils_orig);
                end
            end

            obj.selected_coils_original = sel_orig;

            % Trim padded k-space to the selected coils; from here on coils are renumbered 1..nSel
            kpad = kpad_full(:,:,:,:,:,sel_orig);
            opts = 1:size(kpad,6);   % local coil index space (1..nSel)
            nSel = numel(opts);

            % ----------- FFT to image space (per coil) -----------
            if obj.TwoDimensional == 1
                obj.compleximage = ifft2c(kpad);
            else
                % centered ifft over first 3 spatial dims, keeps all other dims intact
                for t =1:size(kpad,4)
                    for f=1:size(kpad,5)
                        obj.compleximage(:,:,:,t,f,:) = ifft3c(kpad(:,:,:,t,f,:));
                    end
                end
            end


            % ----------- Magnitude k-space (for QC / debugging) -----------
            % Keep in ORIGINAL coil index space for consistency with saved raw k-space
            obj.magkspace = rssq(obj.complexkspace(:,:,:,:,:,sel_orig), 6);

            % ----------- Receiver combination mode (from GUI dropdown) -----------
            combineMode = "ACS/Walsh"; % default
            if isprop(obj,'combination') && ~isempty(obj.combination)
                combineMode = string(obj.combination);
            end

            % ----------- Bias correction mode (from GUI dropdown) -----------
            biasMode = "Off"; % default
            if isprop(obj,'biascorrect') && ~isempty(obj.biascorrect)
                biasMode = string(obj.biascorrect);
            end

            % Initialise outputs
            obj.complexcombined = [];
            obj.magimage = [];
            obj.phaseimage = [];

            % ----------- Combine -----------
            if nSel <= 1
                % Single-coil path: selected coil becomes coil 1 after trimming
                obj.complexcombined = obj.compleximage(:,:,:,:,:,1);
                obj.magimage = abs(obj.complexcombined);

            else
                switch combineMode
                    case "RSSQ"
                        % Magnitude-only
                        obj.complexcombined = [];
                        obj.magimage = rssq(obj.compleximage, 6);

                    case "Adaptive"
                        % Smoothed self-calibrated adaptive complex combine
                        obj.complexcombined = combine_adaptive_smoothed(obj.compleximage, opts, obj);
                        obj.magimage = abs(obj.complexcombined);
                        obj.compleximage = obj.complexcombined;
                    case "ACS/Walsh"
                        % ACS/Walsh sensitivity maps + Roemer combine
                        maps = estimate_maps_from_acs_walsh(kpad, opts, obj);
                        obj.complexcombined = combine_with_maps(obj.compleximage, maps, opts, obj);
                        obj.magimage = abs(obj.complexcombined);
                        obj.compleximage = obj.complexcombined;

                    otherwise
                        
                        obj.complexcombined = [];
                        obj.magimage = rssq(obj.compleximage, 6);
                end

                % ----------- Optional bias-field correction -----------
                switch biasMode
                    case "Off"
                        % no-op

                    case "ACS"
                        [bias, mask3] = bias_field_from_acs(kpad, opts, obj);
                        if exist('apply_bias_field_safe','file') == 2
                            obj.magimage = apply_bias_field_safe(obj.magimage, bias, mask3, obj);
                        else
                            % obj.magimage = apply_bias_field(obj.magimage, bias, mask3, obj);
                            [obj.magimage] = map_guided_bias_correct(obj.magimage, maps, obj);
                        end

                    otherwise
                        % no-op
                end
            end

            % Safety: ensure magimage exists
            if isempty(obj.magimage)
                if ~isempty(obj.complexcombined)
                    obj.magimage = abs(obj.complexcombined);
                else
                    obj.magimage = rssq(obj.compleximage, 6);
                end
            end

            % ----------- Denoise / post-filter -----------
         %  obj.magimage = ffc_mri_filter(obj.magimage, obj.denoise_filter, obj.denoise_params);

            % ----------- Phase image (only meaningful for complex combined) -----------
            if ~isempty(obj.complexcombined)
                obj.phaseimage = angle(obj.complexcombined);
            else
                obj.complexcombined = sum(1.*obj.compleximage,6);
                obj.phaseimage = angle(obj.complexcombined);
            end
        end


        function obj = MapRelaxation(obj)
            %              r1t1mapping;
            dim = size(obj.magimage,1);
            fields = obj.fieldpoints; %list of evolution times in seconds
            times = obj.timepoints;
            n_fields = obj.n_fieldpoints;
            clear obj.T1Maps
            clear obj.R1Maps
            obj.T1Maps = zeros(dim,dim,obj.slices,n_fields(1));
            obj.R1Maps = zeros(dim,dim,obj.slices,n_fields(1));

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
            % obj.T1Maps = T1;
            % obj.R1Maps = R1;
            if obj.checkfit
                for s=1:obj.slices
                    for n=1:n_fields
                        t1map = multipointT1map(squeeze(imagestobeprocessed(:,:,s,:,n)),times(n,:),1,t1mask);
                        T1Maps(:,:,s,n)=t1map(:,:,1,1);
                    end
                end
            else
                for s=1:obj.slices
                    parfor n=1:n_fields
                        t1map = multipointT1map(squeeze(imagestobeprocessed(:,:,s,:,n)),times(n,:),0,t1mask);
                        T1Maps(:,:,s,n)=t1map(:,:,1,1);
                    end
                end
            end


            obj.T1Maps = T1Maps;
            obj.R1Maps = 1000./obj.T1Maps;

        end %image based T1 analysis methods are contained here

 
function obj = T1dispersion(obj, slice)
%T1DISPERSION Complex field-cycling T1 dispersion fit.
%
% Fits the complex ROI signal directly using:
%
%   z(t) = exp(1i*(phi0 + phi1*t)) * (A*exp(-R1*t) + B)
%
% where:
%
%   A    <= 0
%   B    >= 0
%   R1   >= 0
%   phi0 = constant phase offset
%   phi1 = linear phase drift with evolution time, rad/s
%
% Field-cycling interpretation:
%
%   B = M0(B0e)
%   A = -M0(200 mT) - M0(B0e)
%
% This version:
%   - uses obj.compleximage
%   - fits complex data directly
%   - allows linear phase drift across evolution time
%   - tries automatic 180-degree polarity flips for poor initial fits
%   - supports manual polarity flipping in checkfit mode
%   - uses a fast fitter first and a broader fallback only if needed

    % ==========================================================
    % User-adjustable settings
    % ==========================================================
    R1_UPPER_LIMIT = 50;             % s^-1. 50 corresponds to minimum T1 of 20 ms.
    MAX_PHASE_SLOPE = 200;           % rad/s. Increase if low-field phase drift is larger.

    USE_FAST_FITTER_FIRST = true;
    USE_EXPENSIVE_FALLBACK = true;
    FALLBACK_RSQ_THRESHOLD = 0.60;
    USE_PARFOR_STARTLIST = true;

    ENABLE_AUTOMATIC_POLARITY_SEARCH = true;
    AUTO_POLARITY_RSQ_TRIGGER = 0.85;
    AUTO_POLARITY_IMAG_TRIGGER = 0.45;

    MAX_AUTO_POLARITY_FLIPS = 2;
    MAX_AUTO_POLARITY_CANDIDATES = 5;

    AUTO_POLARITY_ACCEPT_RMSE_RATIO = 0.75;
    AUTO_POLARITY_ACCEPT_RSQ_GAIN = 0.15;
    AUTO_POLARITY_MIN_RSQ_ACCEPT = 0.70;

    ENABLE_MANUAL_POLARITY_CORRECTION = true;

    % Keep this off by default. Automatic polarity search is usually more useful.
    ENABLE_AUTOMATIC_TRIMMING = false;

    PRINT_FIT_TIMING = false;

    % ==========================================================
    % Basic setup
    % ==========================================================
    fields = double(obj.fieldpoints(:));      % expected mT
    times  = double(obj.timepoints);          % expected ms

    if isscalar(obj.n_fieldpoints)
        n_fields = double(obj.n_fieldpoints);
    else
        n_fields = double(obj.n_fieldpoints(1));
    end

    n_fields = min(n_fields, numel(fields));

    nParams = 5;                              % [A B R1 phi0 phi1]

    R1dispersion = nan(1, n_fields);
    R1error      = nan(1, n_fields);
    T1dispersion = nan(1, n_fields);
    T1error      = nan(1, n_fields);

    fit_params = nan(n_fields, nParams);      % [A B R1 phi0 phi1]
    fit_sse    = nan(1, n_fields);
    fit_rsq    = nan(1, n_fields);
    fit_rmse   = nan(1, n_fields);
    fit_dfe    = nan(1, n_fields);
    fit_qc     = strings(1, n_fields);
    fit_qc(:)  = "not fitted";

    manual_polarity_flips = cell(1, n_fields);
    auto_polarity_flips   = cell(1, n_fields);
    auto_points_removed   = zeros(1, n_fields);

    last_good_p = [];

    gammaH_Hz_per_T = 42.57747892e6;

    fieldsT  = fields ./ 1e3;
    fieldMHz = fieldsT .* gammaH_Hz_per_T ./ 1e6;

    % ==========================================================
    % Slice-specific mask
    % ==========================================================
    if ndims(obj.mask) == 3
        mask2d = obj.mask(:,:,slice);
    else
        mask2d = obj.mask;
    end

    mask2d = logical(mask2d);

    if ~any(mask2d(:))
        warning('T1dispersion:EmptyMask', ...
            'Mask is empty for slice %d. No fits performed.', slice);
        return
    end

    % ==========================================================
    % Check complex image dimensions
    % Expected: [X Y slice time field]
    % ==========================================================
    ndc = ndims(obj.compleximage);
    csz = size(obj.compleximage);

    if ndc == 6
        error(['T1dispersion: obj.compleximage appears to still have a coil ', ...
               'dimension [X Y slice time field coil]. Perform complex coil ', ...
               'combination before complex T1 fitting.']);
    elseif ndc ~= 5
        error('T1dispersion:UnexpectedComplexImageDims', ...
            'obj.compleximage has unexpected dimensions: %s', mat2str(csz));
    end

    % ==========================================================
    % lsqcurvefit setup
    % ==========================================================
    lsqOpts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'Algorithm', 'trust-region-reflective', ...
        'FunctionTolerance', 1e-7, ...
        'StepTolerance', 1e-7, ...
        'OptimalityTolerance', 1e-7, ...
        'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 1000);

    % Parameter order: [A B R1 phi0 phi1]
    lb = [-Inf, 0, 0,              -pi, -MAX_PHASE_SLOPE];
    ub = [0,    Inf, R1_UPPER_LIMIT, pi,  MAX_PHASE_SLOPE];

    % ==========================================================
    % Main field loop
    % ==========================================================
    for n = 1:n_fields

        if PRINT_FIT_TIMING
            tField = tic;
        end

        % ------------------------------------------------------
        % Extract times for this field
        % ------------------------------------------------------
        if size(times,1) >= n
            invtimes_ms = double(times(n,:)).';
        else
            invtimes_ms = double(times(:));
        end

        invtimes_ms = invtimes_ms(:);

        % ------------------------------------------------------
        % Extract complex image series [X Y T]
        % ------------------------------------------------------
        img_c = reshape( ...
            obj.compleximage(:,:,slice,:,n), ...
            size(obj.compleximage,1), ...
            size(obj.compleximage,2), ...
            size(obj.compleximage,4));

        nT_img = size(img_c, 3);

        if numel(invtimes_ms) ~= nT_img
            nT = min(numel(invtimes_ms), nT_img);

            warning('T1dispersion:TimeImageMismatch', ...
                'Field %d: time vector length does not match image time dimension. Using first %d points.', ...
                n, nT);

            invtimes_ms = invtimes_ms(1:nT);
            img_c = img_c(:,:,1:nT);
        end

        valid_time = isfinite(invtimes_ms);

        invtimes_ms = invtimes_ms(valid_time);
        img_c = img_c(:,:,valid_time);

        if numel(invtimes_ms) < 4
            warning('T1dispersion:TooFewTimepoints', ...
                'Field %d has fewer than 4 valid timepoints. Skipping.', n);
            fit_qc(n) = "too few points";
            continue
        end

        % ------------------------------------------------------
        % Complex ROI signal
        % ------------------------------------------------------
        z_roi = complexRoiSignal(img_c, mask2d);

        valid = isfinite(invtimes_ms) & ...
                isfinite(real(z_roi)) & ...
                isfinite(imag(z_roi));

        x_ms = invtimes_ms(valid);
        z    = z_roi(valid);

        if numel(z) < 4
            warning('T1dispersion:TooFewValidSignalPoints', ...
                'Field %d has fewer than 4 valid complex signal points. Skipping.', n);
            fit_qc(n) = "too few valid points";
            continue
        end

        [x_ms, ord] = sort(x_ms);
        z = z(ord);

        % Previous R1 as a useful starting point
        if n > 1 && isfinite(R1dispersion(n-1)) && R1dispersion(n-1) > 0
            R1start_centre = R1dispersion(n-1);
        else
            R1start_centre = 1; % s^-1
        end

        % ------------------------------------------------------
        % Initial complex fit
        % ------------------------------------------------------
        [p, ci, gof, ok, xs, z_used, residual_complex] = fitComplexFcModelAuto( ...
            x_ms, z, R1start_centre, lb, ub, lsqOpts, last_good_p);

        if ~ok
            warning('T1dispersion:FitFailed', ...
                'Complex fit failed for field %d.', n);
            fit_qc(n) = "fit failed";
            continue
        end

        % ------------------------------------------------------
        % Automatic polarity search
        % ------------------------------------------------------
        auto_flip_idx = [];

        if ENABLE_AUTOMATIC_POLARITY_SEARCH && shouldTryAutoPolaritySearch(gof, p, xs, z_used)

            [p_auto, ci_auto, gof_auto, ok_auto, xs_auto, z_auto, residual_auto, flip_idx_auto] = ...
                automaticPolarityFlipSearch( ...
                    p, ci, gof, xs, z_used, residual_complex, ...
                    R1start_centre, lb, ub, lsqOpts);

            if ok_auto && ~isempty(flip_idx_auto)

                rmseImproved = isfinite(gof_auto.rmse) && isfinite(gof.rmse) && ...
                    gof_auto.rmse < AUTO_POLARITY_ACCEPT_RMSE_RATIO .* gof.rmse;

                rsqImproved = isfinite(gof_auto.rsquare) && isfinite(gof.rsquare) && ...
                    gof_auto.rsquare > gof.rsquare + AUTO_POLARITY_ACCEPT_RSQ_GAIN && ...
                    gof_auto.rsquare > AUTO_POLARITY_MIN_RSQ_ACCEPT;

                if rmseImproved || rsqImproved

                    p = p_auto;
                    ci = ci_auto;
                    gof = gof_auto;
                    xs = xs_auto;
                    z_used = z_auto;
                    residual_complex = residual_auto;

                    auto_flip_idx = flip_idx_auto;
                    auto_polarity_flips{n} = auto_flip_idx;

                    warning('T1dispersion:AutoPolarityFlip', ...
                        'Field %d: automatically flipped point(s) %s. New R^2 = %.4f.', ...
                        n, mat2str(auto_flip_idx), gof.rsquare);
                end
            end
        end

        % ------------------------------------------------------
        % Optional automatic trimming: disabled by default
        % ------------------------------------------------------
        if ENABLE_AUTOMATIC_TRIMMING && ok && numel(z_used) >= 6 && gof.rsquare < 0.70
            [p2, ci2, gof2, ok2, xs2, z2, residual2, removed_idx] = ...
                fitComplexFcModelWithSingleTrim( ...
                    xs .* 1000, z_used, R1start_centre, lb, ub, lsqOpts, p);

            if ok2 && isfinite(gof2.rmse) && isfinite(gof.rmse) && gof2.rmse < 0.75 .* gof.rmse
                p = p2;
                ci = ci2;
                gof = gof2;
                xs = xs2;
                z_used = z2;
                residual_complex = residual2;
                auto_points_removed(n) = numel(removed_idx);

                warning('T1dispersion:AutoPointRemoved', ...
                    'Field %d: automatically removed point(s): %s', ...
                    n, mat2str(removed_idx));
            end
        end

        % ------------------------------------------------------
        % Optional manual polarity correction during checkfit
        % ------------------------------------------------------
        manual_flip_idx = [];

        if ENABLE_MANUAL_POLARITY_CORRECTION && hasProp(obj, 'checkfit') && obj.checkfit == 1

            [p, ci, gof, xs, z_used, residual_complex, manual_flip_idx] = ...
                manualPolarityCorrectionLoop( ...
                    n, p, ci, gof, xs, z_used, residual_complex, ...
                    R1start_centre, lb, ub, lsqOpts, auto_flip_idx);

            manual_polarity_flips{n} = manual_flip_idx;

        elseif hasProp(obj, 'checkfit') && obj.checkfit == 1

            showStaticFitFigure(n, p, gof, xs, z_used, auto_flip_idx);

        end

        % ------------------------------------------------------
        % Store result
        % ------------------------------------------------------
        A_fit    = p(1); %#ok<NASGU>
        B_fit    = p(2); %#ok<NASGU>
        R1_fit   = p(3);
        phi0_fit = p(4); %#ok<NASGU>
        phi1_fit = p(5); %#ok<NASGU>

        R1dispersion(n) = R1_fit;
        fit_params(n,:) = p;

        fit_sse(n)  = gof.sse;
        fit_rsq(n)  = gof.rsquare;
        fit_rmse(n) = gof.rmse;
        fit_dfe(n)  = gof.dfe;

        if all(isfinite(ci(:))) && size(ci,1) == 2 && size(ci,2) >= 3
            R1error(n) = diff(ci(:,3)) ./ 2;
        else
            R1error(n) = nan;
        end

        if gof.rsquare < 0.70
            fit_qc(n) = "poor";
        elseif gof.rsquare < 0.90
            fit_qc(n) = "questionable";
        else
            fit_qc(n) = "good";
        end

        if ~isempty(auto_flip_idx)
            fit_qc(n) = fit_qc(n) + " / auto-flipped";
        end

        if ~isempty(manual_flip_idx)
            fit_qc(n) = fit_qc(n) + " / manually corrected";
        end

        if isfinite(R1_fit) && R1_fit > 0
            last_good_p = p;
        end

        if PRINT_FIT_TIMING
            fprintf('Field %d fit time: %.3f s, R1 = %.4g s^-1, R2 = %.4f\n', ...
                n, toc(tField), R1_fit, gof.rsquare);
        end

    end

    % ==========================================================
    % Convert R1 to T1 with proper uncertainty propagation
    % ==========================================================
    valid_R1 = isfinite(R1dispersion) & R1dispersion > 0;

    T1dispersion(valid_R1) = 1000 ./ R1dispersion(valid_R1);

    % T1 = 1000/R1
    % dT1 = 1000*dR1/R1^2
    T1error(valid_R1) = 1000 .* R1error(valid_R1) ./ ...
        (R1dispersion(valid_R1).^2);

    % ==========================================================
    % Display and plot
    % ==========================================================
    if hasProp(obj, 'R1T1') && obj.R1T1 == 1

        disp(['Fields (MHz): ' num2str(fieldMHz(1:n_fields).')])
        disp(['R1 (s^-1): ' num2str(R1dispersion)])
        disp(['R1 error (s^-1): ' num2str(R1error)])
        disp(['Fit QC: ' char(join(fit_qc, ', '))])

        valid_plot = isfinite(fieldMHz(1:n_fields)) & fieldMHz(1:n_fields) > 0 & ...
                     isfinite(R1dispersion(:).') & R1dispersion > 0;

        figure;
        scatter(fieldMHz(valid_plot), R1dispersion(valid_plot), 'o');
        set(gca, 'xscale', 'log', 'yscale', 'log');
        xlabel('Evolution Field (MHz)');
        ylabel('R_1 (s^{-1})');
        title('Complex T1 dispersion fit: R_1 versus evolution field');
        grid on;

    else

        disp(['Fields (T): ' num2str(fieldsT(1:n_fields).')])
        disp(['T1 (ms): ' num2str(T1dispersion)])
        disp(['T1 error (ms): ' num2str(T1error)])
        disp(['R1 (s^-1): ' num2str(R1dispersion)])
        disp(['Fit QC: ' char(join(fit_qc, ', '))])

        valid_plot = isfinite(fieldsT(1:n_fields)) & fieldsT(1:n_fields) > 0 & ...
                     isfinite(T1dispersion(:).') & T1dispersion > 0;

        figure;
        plot(fieldsT(valid_plot), T1dispersion(valid_plot), 'o');
        set(gca, 'xscale', 'log', 'yscale', 'log');
        xlabel('Evolution Field (T)');
        ylabel('T_1 (ms)');
        title('Complex T1 dispersion fit: T_1 versus evolution field');
        grid on;

    end

    % ==========================================================
    % Store results where possible
    % ==========================================================
    obj = safeSet(obj, 'dispersioncurve', R1dispersion);

    obj = safeSet(obj, 'R1dispersion', R1dispersion);
    obj = safeSet(obj, 'R1error', R1error);
    obj = safeSet(obj, 'T1dispersion', T1dispersion);
    obj = safeSet(obj, 'T1error', T1error);

    obj = safeSet(obj, 'fit_params', fit_params);
    obj = safeSet(obj, 'fit_sse', fit_sse);
    obj = safeSet(obj, 'fit_rsq', fit_rsq);
    obj = safeSet(obj, 'fit_rmse', fit_rmse);
    obj = safeSet(obj, 'fit_dfe', fit_dfe);
    obj = safeSet(obj, 'fit_qc', fit_qc);

    obj = safeSet(obj, 'manual_polarity_flips', manual_polarity_flips);
    obj = safeSet(obj, 'auto_polarity_flips', auto_polarity_flips);
    obj = safeSet(obj, 'auto_points_removed', auto_points_removed);

    % ==========================================================
    % Nested helper: complex ROI signal
    % ==========================================================
    function z_roi = complexRoiSignal(img_c, mask2d)

        nT_local = size(img_c, 3);
        z_roi = nan(nT_local, 1);

        for tt = 1:nT_local
            tmp = img_c(:,:,tt);
            vals = tmp(mask2d);

            finiteVals = isfinite(real(vals)) & isfinite(imag(vals));
            vals = vals(finiteVals);

            if isempty(vals)
                z_roi(tt) = nan;
            else
                z_roi(tt) = mean(vals);
            end
        end

    end

    % ==========================================================
    % Nested helper: automatic fitter
    % ==========================================================
    function [p, ci, gof, ok, xs, z_used, residual_complex] = fitComplexFcModelAuto( ...
            x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint)

        if USE_FAST_FITTER_FIRST
            [p, ci, gof, ok, xs, z_used, residual_complex] = fitComplexFcModelFast( ...
                x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint);

            needsFallback = ~ok || ~isfinite(gof.rsquare) || gof.rsquare < FALLBACK_RSQ_THRESHOLD;

            if USE_EXPENSIVE_FALLBACK && needsFallback
                [p2, ci2, gof2, ok2, xs2, z2, residual2] = fitComplexFcModelBroad( ...
                    x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint);

                if ok2 && (~ok || gof2.sse < gof.sse)
                    p = p2;
                    ci = ci2;
                    gof = gof2;
                    ok = ok2;
                    xs = xs2;
                    z_used = z2;
                    residual_complex = residual2;
                end
            end
        else
            [p, ci, gof, ok, xs, z_used, residual_complex] = fitComplexFcModelBroad( ...
                x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint);
        end

    end

    % ==========================================================
    % Nested helper: fast complex fitter
    % ==========================================================
    function [bestP, bestCi, bestGof, ok, xs, z_used, residual_complex] = fitComplexFcModelFast( ...
            x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint)

        ok = false;

        bestP = nan(1, numel(lb));
        bestCi = nan(2, numel(lb));
        bestGof = emptyGof(Inf);

        x_s = double(x_ms(:)) ./ 1000;
        z   = complex(z(:));

        [xs, ord] = sort(x_s);
        z_used = z(ord);

        ydata = [real(z_used); imag(z_used)];
        residual_complex = nan(size(z_used));

        % Initial phase from early high-SNR point.
        n_early = min(3, numel(z_used));
        early_idx = 1:n_early;

        [~, idx_local] = max(abs(z_used(early_idx)));
        ref_idx = early_idx(idx_local);

        phi0 = wrapPiLocal(angle(z_used(ref_idx)) - pi);

        phi0Starts = wrapPiLocal(phi0 + [-pi/12, 0, pi/12]);

        if ~isempty(p_hint) && numel(p_hint) >= 4 && isfinite(p_hint(4))
            phi0Starts = [p_hint(4), phi0Starts];
        end

        phi0Starts = unique(round(wrapPiLocal(phi0Starts), 10), 'stable');

        % Include small phase slopes. Optimiser can move from these.
        phi1Starts = [0, -10, 10];

        if ~isempty(p_hint) && numel(p_hint) >= 5 && isfinite(p_hint(5))
            phi1Starts = [p_hint(5), phi1Starts];
        end

        phi1Starts = phi1Starts(phi1Starts >= lb(5) & phi1Starts <= ub(5));
        phi1Starts = unique(round(phi1Starts, 10), 'stable');

        R1start_centre = max(R1start_centre, 0.001);

        if ~isempty(p_hint) && numel(p_hint) >= 3 && isfinite(p_hint(3)) && p_hint(3) > 0
            R1base = p_hint(3);
        else
            R1base = R1start_centre;
        end

        R1starts = [ ...
            R1base, ...
            R1base/2, ...
            R1base*2, ...
            R1start_centre, ...
            0.1, 0.5, 1, 2, 5, 10];

        R1starts = R1starts(isfinite(R1starts) & R1starts > 0 & R1starts <= ub(3));
        R1starts = unique(round(R1starts, 10), 'stable');

        % Build start list
        p0_list = [];

        p_hint = normaliseParamVector(p_hint);

        if ~isempty(p_hint) && all(isfinite(p_hint))
            p0 = max(p_hint, lb);
            p0 = min(p0, ub);
            p0_list = [p0_list; p0]; %#ok<AGROW>
        end

        for pp = 1:numel(phi0Starts)
            for qq = 1:numel(phi1Starts)

                phi0_start = phi0Starts(pp);
                phi1_start = phi1Starts(qq);

                phase_start = phi0_start + phi1_start .* xs;
                yrot = real(z_used .* exp(-1i .* phase_start));

                amp_scale = max(abs(yrot));
                if ~isfinite(amp_scale) || amp_scale == 0
                    amp_scale = max(abs(z_used));
                end
                if ~isfinite(amp_scale) || amp_scale == 0
                    amp_scale = 1;
                end

                smallNeg = -max(1e-9 .* amp_scale, 1e-12);

                for rr = 1:numel(R1starts)

                    R1_start = R1starts(rr);

                    u = exp(-xs .* R1_start);
                    X = [u, ones(size(u))];

                    beta = X \ yrot;

                    A0 = beta(1);
                    B0 = beta(2);

                    A0 = min(A0, smallNeg);
                    B0 = max(B0, 0);

                    p0 = [A0, B0, R1_start, phi0_start, phi1_start];

                    p0 = max(p0, lb);
                    p0 = min(p0, ub);

                    p0_list = [p0_list; p0]; %#ok<AGROW>
                end
            end
        end

        if ~isempty(p0_list)
            p0_list = unique(round(p0_list, 8), 'rows', 'stable');
        end

        [bestP, bestCi, bestGof, ok, residual_complex] = runStartList( ...
            p0_list, xs, z_used, ydata, lb, ub, lsqOpts, bestP, bestCi, bestGof, ok);

    end

    % ==========================================================
    % Nested helper: broader fallback complex fitter
    % ==========================================================
    function [bestP, bestCi, bestGof, ok, xs, z_used, residual_complex] = fitComplexFcModelBroad( ...
            x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint)

        ok = false;

        bestP = nan(1, numel(lb));
        bestCi = nan(2, numel(lb));
        bestGof = emptyGof(Inf);

        x_s = double(x_ms(:)) ./ 1000;
        z   = complex(z(:));

        [xs, ord] = sort(x_s);
        z_used = z(ord);

        ydata = [real(z_used); imag(z_used)];
        residual_complex = nan(size(z_used));

        n_early = min(3, numel(z_used));
        early_idx = 1:n_early;

        [~, idx_local] = max(abs(z_used(early_idx)));
        ref_idx = early_idx(idx_local);

        phi0 = wrapPiLocal(angle(z_used(ref_idx)) - pi);

        phi0Starts = wrapPiLocal([ ...
            phi0, ...
            phi0 - pi/2, phi0 - pi/4, phi0 - pi/8, ...
            phi0 + pi/8, phi0 + pi/4, phi0 + pi/2, ...
            angle(mean(z_used))]);

        if ~isempty(p_hint) && numel(p_hint) >= 4 && isfinite(p_hint(4))
            phi0Starts = [p_hint(4), phi0Starts];
        end

        phi0Starts = unique(round(wrapPiLocal(phi0Starts), 10), 'stable');

        phi1Starts = [0, -5, 5, -20, 20, -50, 50];

        if ~isempty(p_hint) && numel(p_hint) >= 5 && isfinite(p_hint(5))
            phi1Starts = [p_hint(5), phi1Starts];
        end

        phi1Starts = phi1Starts(phi1Starts >= lb(5) & phi1Starts <= ub(5));
        phi1Starts = unique(round(phi1Starts, 10), 'stable');

        R1start_centre = max(R1start_centre, 0.001);

        R1starts = [ ...
            R1start_centre/5, ...
            R1start_centre/2, ...
            R1start_centre, ...
            R1start_centre*2, ...
            R1start_centre*5, ...
            0.005, 0.01, 0.02, 0.05, ...
            0.1, 0.2, 0.5, 1, 2, 5, 10, 20];

        if ~isempty(p_hint) && numel(p_hint) >= 3 && isfinite(p_hint(3)) && p_hint(3) > 0
            R1starts = [p_hint(3), R1starts];
        end

        R1starts = R1starts(isfinite(R1starts) & R1starts > 0 & R1starts <= ub(3));
        R1starts = unique(round(R1starts, 10), 'stable');

        p0_list = [];

        p_hint = normaliseParamVector(p_hint);

        if ~isempty(p_hint) && all(isfinite(p_hint))
            p0 = max(p_hint, lb);
            p0 = min(p0, ub);
            p0_list = [p0_list; p0]; %#ok<AGROW>
        end

        for pp = 1:numel(phi0Starts)
            for qq = 1:numel(phi1Starts)

                phi0_start = phi0Starts(pp);
                phi1_start = phi1Starts(qq);

                phase_start = phi0_start + phi1_start .* xs;
                yrot = real(z_used .* exp(-1i .* phase_start));

                amp_scale = max(abs(yrot));
                if ~isfinite(amp_scale) || amp_scale == 0
                    amp_scale = max(abs(z_used));
                end
                if ~isfinite(amp_scale) || amp_scale == 0
                    amp_scale = 1;
                end

                smallNeg = -max(1e-9 .* amp_scale, 1e-12);

                y_first = mean(yrot(1:min(2,end)));
                y_last  = mean(yrot(max(1,end-1):end));

                B_guess_1 = max(y_last, 0);
                B_guess_2 = max(max(yrot), 0);
                B_guess_3 = max(median(yrot), 0);

                Bstarts = unique([B_guess_1, B_guess_2, B_guess_3, 0], 'stable');

                A_guess_1 = min(y_first - B_guess_1, smallNeg);
                A_guess_2 = min(min(yrot) - B_guess_2, smallNeg);
                A_guess_3 = min(-amp_scale, smallNeg);
                A_guess_4 = min(-2 .* amp_scale, smallNeg);

                Astarts = unique([A_guess_1, A_guess_2, A_guess_3, A_guess_4], 'stable');

                for aa = 1:numel(Astarts)
                    for bb = 1:numel(Bstarts)
                        for rr = 1:numel(R1starts)

                            p0 = [Astarts(aa), Bstarts(bb), R1starts(rr), phi0_start, phi1_start];

                            p0 = max(p0, lb);
                            p0 = min(p0, ub);

                            p0_list = [p0_list; p0]; %#ok<AGROW>
                        end
                    end
                end
            end
        end

        if ~isempty(p0_list)
            p0_list = unique(round(p0_list, 8), 'rows', 'stable');
        end

        [bestP, bestCi, bestGof, ok, residual_complex] = runStartList( ...
            p0_list, xs, z_used, ydata, lb, ub, lsqOpts, bestP, bestCi, bestGof, ok);

    end

    % ==========================================================
    % Nested helper: run lsqcurvefit from start list
    % ==========================================================
   function [bestP, bestCi, bestGof, ok, residual_complex] = runStartList( ...
        p0_list, xs, z_used, ydata, lb, ub, lsqOpts, bestP, bestCi, bestGof, ok)

    residual_complex = nan(size(z_used));

    if isempty(p0_list)
        return
    end

    nStarts = size(p0_list, 1);
    nParams = numel(lb);
    nZ = numel(z_used);

    % Only use parfor when there are enough starts to justify the overhead.
    usePar = USE_PARFOR_STARTLIST && nStarts >= 8;

    if usePar

        p_all        = nan(nStarts, nParams);
        sse_all      = inf(nStarts, 1);
        rsq_all      = nan(nStarts, 1);
        adjrsq_all   = nan(nStarts, 1);
        dfe_all      = nan(nStarts, 1);
        rmse_all     = nan(nStarts, 1);
        exitflag_all = nan(nStarts, 1);
        residual_all = nan(numel(ydata), nStarts);

        % Important:
        % Do not call nested functions inside this parfor.
        parfor ii = 1:nStarts

            p0 = p0_list(ii,:);

            local_p        = nan(1, nParams);
            local_sse      = inf;
            local_rsq      = nan;
            local_adjrsq   = nan;
            local_dfe      = nan;
            local_rmse     = nan;
            local_exitflag = nan;
            local_residual = nan(numel(ydata), 1);

            try
                % Self-contained model. Do NOT call complexModelVector here.
                localModelFun = @(p, x) [ ...
                    real(exp(1i .* (p(4) + p(5) .* x(:))) .* ...
                         (p(1) .* exp(-x(:) .* p(3)) + p(2))); ...
                    imag(exp(1i .* (p(4) + p(5) .* x(:))) .* ...
                         (p(1) .* exp(-x(:) .* p(3)) + p(2)))];

                [pfit, resnorm, residual, exitflag] = ...
                    lsqcurvefit(localModelFun, p0, xs, ydata, lb, ub, lsqOpts);

                if isfinite(resnorm)

                    local_p = pfit(:).';
                    local_sse = resnorm;
                    local_exitflag = exitflag;
                    local_residual = residual(:);

                    nObs = numel(ydata);
                    nPar = numel(local_p);
                    dfe = max(nObs - nPar, 0);

                    sst = sum((ydata - mean(ydata)).^2);

                    if sst > 0
                        rsq = 1 - resnorm ./ sst;
                    else
                        rsq = NaN;
                    end

                    if dfe > 0 && isfinite(rsq)
                        adjrsq = 1 - (1 - rsq) .* ((nObs - 1) ./ dfe);
                        rmse = sqrt(resnorm ./ dfe);
                    else
                        adjrsq = NaN;
                        rmse = NaN;
                    end

                    local_rsq = rsq;
                    local_adjrsq = adjrsq;
                    local_dfe = dfe;
                    local_rmse = rmse;

                end

            catch
                % Leave local result as failed.
            end

            p_all(ii,:)        = local_p;
            sse_all(ii)        = local_sse;
            rsq_all(ii)        = local_rsq;
            adjrsq_all(ii)     = local_adjrsq;
            dfe_all(ii)        = local_dfe;
            rmse_all(ii)       = local_rmse;
            exitflag_all(ii)   = local_exitflag;
            residual_all(:,ii) = local_residual;

        end

        [bestSse, bestIdx] = min(sse_all);

        if isfinite(bestSse)

            bestP = p_all(bestIdx,:);

            bestGof = struct( ...
                'sse',        sse_all(bestIdx), ...
                'rsquare',    rsq_all(bestIdx), ...
                'adjrsquare', adjrsq_all(bestIdx), ...
                'dfe',        dfe_all(bestIdx), ...
                'rmse',       rmse_all(bestIdx), ...
                'exitflag',   exitflag_all(bestIdx));

            residual = residual_all(:,bestIdx);
            residual_complex = residual(1:nZ) + 1i .* residual(nZ+1:end);

            ok = true;

            % Re-run only the best solution once, in serial, to get the Jacobian
            % for confidence intervals. This avoids calling nested helper
            % functions inside parfor.
            try
                serialModelFun = @(p, x) complexModelVector(p, x);

                [pfit2, resnorm2, residual2, exitflag2, ~, ~, jacobian2] = ...
                    lsqcurvefit(serialModelFun, bestP, xs, ydata, lb, ub, lsqOpts);

                if isfinite(resnorm2) && resnorm2 <= bestGof.sse .* 1.001

                    bestP = pfit2(:).';

                    nObs = numel(ydata);
                    nPar = numel(bestP);
                    dfe = max(nObs - nPar, 0);

                    sst = sum((ydata - mean(ydata)).^2);

                    if sst > 0
                        rsq = 1 - resnorm2 ./ sst;
                    else
                        rsq = NaN;
                    end

                    if dfe > 0 && isfinite(rsq)
                        adjrsq = 1 - (1 - rsq) .* ((nObs - 1) ./ dfe);
                        rmse = sqrt(resnorm2 ./ dfe);
                    else
                        adjrsq = NaN;
                        rmse = NaN;
                    end

                    bestGof = struct( ...
                        'sse',        resnorm2, ...
                        'rsquare',    rsq, ...
                        'adjrsquare', adjrsq, ...
                        'dfe',        dfe, ...
                        'rmse',       rmse, ...
                        'exitflag',   exitflag2);

                    bestCi = approximateParameterCi(bestP, resnorm2, jacobian2, dfe);

                    residual_complex = residual2(1:nZ) + 1i .* residual2(nZ+1:end);

                else
                    bestCi = nan(2, nParams);
                end

            catch
                bestCi = nan(2, nParams);
            end

        end

    else

        % Serial version can still call nested functions.
        modelFun = @(p, x) complexModelVector(p, x);

        for ii = 1:nStarts

            p0 = p0_list(ii,:);

            try
                [pfit, resnorm, residual, exitflag, ~, ~, jacobian] = ...
                    lsqcurvefit(modelFun, p0, xs, ydata, lb, ub, lsqOpts);

                if isfinite(resnorm) && resnorm < bestGof.sse

                    bestP = pfit(:).';
                    ok = true;

                    nObs = numel(ydata);
                    nPar = numel(bestP);
                    dfe = max(nObs - nPar, 0);

                    sst = sum((ydata - mean(ydata)).^2);

                    if sst > 0
                        rsq = 1 - resnorm ./ sst;
                    else
                        rsq = NaN;
                    end

                    if dfe > 0 && isfinite(rsq)
                        adjrsq = 1 - (1 - rsq) .* ((nObs - 1) ./ dfe);
                        rmse = sqrt(resnorm ./ dfe);
                    else
                        adjrsq = NaN;
                        rmse = NaN;
                    end

                    bestGof = struct( ...
                        'sse',        resnorm, ...
                        'rsquare',    rsq, ...
                        'adjrsquare', adjrsq, ...
                        'dfe',        dfe, ...
                        'rmse',       rmse, ...
                        'exitflag',   exitflag);

                    bestCi = approximateParameterCi(bestP, resnorm, jacobian, dfe);

                    residual_complex = residual(1:nZ) + 1i .* residual(nZ+1:end);

                end

            catch
                % Try next start point.
            end
        end
    end
end

    % ==========================================================
    % Nested helper: complex model vector
    % ==========================================================
    function yvec = complexModelVector(p, x_s)
        % p = [A B R1 phi0 phi1]

        A    = p(1);
        B    = p(2);
        R1   = p(3);
        phi0 = p(4);
        phi1 = p(5);

        x_s = x_s(:);

        realSignal = A .* exp(-x_s .* R1) + B;
        phase = phi0 + phi1 .* x_s;

        zModel = exp(1i .* phase) .* realSignal;

        yvec = [real(zModel); imag(zModel)];
    end

    % ==========================================================
    % Nested helper: should automatic polarity search run?
    % ==========================================================
    function tf = shouldTryAutoPolaritySearch(gof, p, xs, z_used)

        tf = false;

        if ~isfinite(gof.rsquare) || gof.rsquare < AUTO_POLARITY_RSQ_TRIGGER
            tf = true;
            return
        end

        [imagFrac, ~] = imaginaryFractionAfterFit(p, xs, z_used);

        if isfinite(imagFrac) && imagFrac > AUTO_POLARITY_IMAG_TRIGGER
            tf = true;
        end

    end

    % ==========================================================
    % Nested helper: automatic polarity flip search
    % ==========================================================
    function [bestP, bestCi, bestGof, okBest, bestXs, bestZ, bestResidual, bestFlipIdx] = ...
        automaticPolarityFlipSearch( ...
            p_current, ci_current, gof_current, xs_current, z_current, residual_current, ...
            R1start_centre, lb, ub, lsqOpts)

        bestP = p_current;
        bestCi = ci_current;
        bestGof = gof_current;
        bestXs = xs_current;
        bestZ = z_current;
        bestResidual = residual_current;
        bestFlipIdx = [];
        okBest = true;

        nPts = numel(z_current);

        if nPts < 4
            return
        end

        maxFlips = min(MAX_AUTO_POLARITY_FLIPS, nPts);
        maxCandidates = min(MAX_AUTO_POLARITY_CANDIDATES, nPts);

        candidateIdx = rankPolarityFlipCandidates(p_current, xs_current, z_current, maxCandidates);

        if isempty(candidateIdx)
            return
        end

        baselineScore = autoPolarityScore(gof_current, 0, p_current, xs_current, z_current, ub);

        for nFlip = 1:maxFlips

            if numel(candidateIdx) < nFlip
                continue
            end

            combos = nchoosek(candidateIdx(:).', nFlip);

            for cc = 1:size(combos, 1)

                flipIdx = combos(cc,:);

                z_try = z_current;
                z_try(flipIdx) = -z_try(flipIdx);

                % Use the fast fitter only here to keep the automatic search manageable.
                [p_try, ci_try, gof_try, ok_try, xs_try, z_used_try, residual_try] = ...
                    fitComplexFcModelFast( ...
                        xs_current .* 1000, z_try, ...
                        R1start_centre, lb, ub, lsqOpts, p_current);

                if ~ok_try
                    continue
                end

                tryScore = autoPolarityScore(gof_try, nFlip, p_try, xs_try, z_used_try, ub);
                bestScore = autoPolarityScore(bestGof, numel(bestFlipIdx), bestP, bestXs, bestZ, ub);

                betterThanBest = isfinite(tryScore) && tryScore < bestScore;
                betterThanBaseline = isfinite(tryScore) && tryScore < baselineScore;

                if betterThanBest && betterThanBaseline
                    bestP = p_try;
                    bestCi = ci_try;
                    bestGof = gof_try;
                    bestXs = xs_try;
                    bestZ = z_used_try;
                    bestResidual = residual_try;
                    bestFlipIdx = sort(flipIdx);
                    okBest = true;
                end
            end
        end

    end

    % ==========================================================
    % Nested helper: rank likely polarity-error points
    % ==========================================================
    function candidateIdx = rankPolarityFlipCandidates(p, xs, z_used, maxCandidates)

        candidateIdx = [];

        if isempty(z_used) || isempty(xs) || numel(p) < 5
            return
        end

        A    = p(1);
        B    = p(2);
        R1   = p(3);
        phi0 = p(4);
        phi1 = p(5);

        xs = xs(:);
        z_used = z_used(:);

        phase_fit = phi0 + phi1 .* xs;

        z_rot = z_used .* exp(-1i .* phase_fit);

        y_real = real(z_rot);
        y_imag = imag(z_rot);

        y_fit = A .* exp(-xs .* R1) + B;

        scale = max(abs([y_real(:); y_fit(:)]));
        if ~isfinite(scale) || scale == 0
            scale = max(abs(z_used));
        end
        if ~isfinite(scale) || scale == 0
            scale = 1;
        end

        realResidual = abs(y_real - y_fit) ./ scale;
        imagResidual = abs(y_imag) ./ scale;

        wrongSign = sign(y_real) ~= sign(y_fit) & ...
                    abs(y_real) > 0.15 .* scale & ...
                    abs(y_fit)  > 0.15 .* scale;

        candidateScore = realResidual + 0.75 .* imagResidual + 1.0 .* double(wrongSign);

        nearZero = abs(y_real) < 0.05 .* scale & abs(y_fit) < 0.05 .* scale;
        candidateScore(nearZero) = 0.5 .* candidateScore(nearZero);

        finiteScore = isfinite(candidateScore);

        if ~any(finiteScore)
            return
        end

        candidateScore(~finiteScore) = -Inf;

        [~, order] = sort(candidateScore, 'descend');

        nTake = min(maxCandidates, numel(order));
        candidateIdx = order(1:nTake).';

    end

    % ==========================================================
    % Nested helper: score automatic polarity candidate
    % ==========================================================
    function score = autoPolarityScore(gof, nFlips, p, xs, z_used, ub)

        score = Inf;

        if isempty(gof) || ~isfinite(gof.rmse)
            return
        end

        score = gof.rmse;

        % Penalise flipping more points.
        score = score .* (1 + 0.20 .* nFlips);

        % Penalise fits that hit parameter bounds.
        if numel(p) >= 5

            R1 = p(3);
            phi1 = p(5);

            if isfinite(R1) && R1 > 0.95 .* ub(3)
                score = score .* 1.25;
            end

            if isfinite(phi1) && abs(phi1) > 0.95 .* ub(5)
                score = score .* 1.20;
            end
        end

        [imagFrac, ~] = imaginaryFractionAfterFit(p, xs, z_used);

        if isfinite(imagFrac)
            score = score .* (1 + 0.5 .* imagFrac);
        end

    end

    % ==========================================================
    % Nested helper: imaginary-component diagnostic
    % ==========================================================
    function [imagFrac, scale] = imaginaryFractionAfterFit(p, xs, z_used)

        imagFrac = NaN;
        scale = NaN;

        if isempty(z_used) || isempty(xs) || numel(p) < 5
            return
        end

        phi0 = p(4);
        phi1 = p(5);

        phase_fit = phi0 + phi1 .* xs(:);

        z_rot = z_used(:) .* exp(-1i .* phase_fit);

        y_real = real(z_rot);
        y_imag = imag(z_rot);

        scale = median(abs(y_real), 'omitnan');

        if ~isfinite(scale) || scale == 0
            scale = median(abs(z_used), 'omitnan');
        end

        if ~isfinite(scale) || scale == 0
            return
        end

        imagFrac = median(abs(y_imag), 'omitnan') ./ scale;

    end

    % ==========================================================
    % Nested helper: manual polarity correction loop
    % ==========================================================
    function [p, ci, gof, xs, z_used, residual_complex, manual_flip_idx] = ...
        manualPolarityCorrectionLoop( ...
            field_index, p, ci, gof, xs, z_used, residual_complex, ...
            R1start_centre, lb, ub, lsqOpts, auto_flip_idx)

        manual_flip_idx = [];

        h = figure('Name', sprintf('Manual polarity correction: field %d', field_index));

        while true

            drawManualComplexFitFigure( ...
                h, field_index, xs, z_used, p, gof, manual_flip_idx, auto_flip_idx);

            fprintf('\nField %d manual correction:\n', field_index);
            fprintf('  Left-click a point to flip its polarity and refit.\n');
            fprintf('  Press Enter or right-click to accept and continue.\n');
            fprintf('  Press u to undo the most recent manual flip.\n');

            figure(h);

            try
                [xclick, yclick, button] = ginput(1);
            catch
                break
            end

            if isempty(button)
                break
            end

            if button == 3
                break
            end

            % Keyboard undo
            if button == double('u') || button == double('U')
                if ~isempty(manual_flip_idx)

                    undo_idx = manual_flip_idx(end);
                    z_used(undo_idx) = -z_used(undo_idx);
                    manual_flip_idx(end) = [];

                    R1start_refit = p(3);
                    if ~isfinite(R1start_refit) || R1start_refit <= 0
                        R1start_refit = R1start_centre;
                    end

                    [p2, ci2, gof2, ok2, xs2, z2, residual2] = fitComplexFcModelAuto( ...
                        xs .* 1000, z_used, R1start_refit, lb, ub, lsqOpts, p);

                    if ok2
                        p = p2;
                        ci = ci2;
                        gof = gof2;
                        xs = xs2;
                        z_used = z2;
                        residual_complex = residual2;
                    else
                        z_used(undo_idx) = -z_used(undo_idx);
                        manual_flip_idx(end+1) = undo_idx;

                        warndlg('Undo caused the refit to fail, so it was reverted.', ...
                            'Refit failed');
                    end
                end

                continue
            end

            % Only left-click flips a point
            if button ~= 1
                continue
            end

            idx = nearestRealFitPoint(xs, z_used, p, xclick, yclick);

            if isempty(idx) || ~isfinite(idx)
                continue
            end

            msg = sprintf(['Flip polarity of point %d?\n\n', ...
                           'Evolution time = %.4g s\n', ...
                           'This will multiply the complex ROI point by -1 and refit.'], ...
                           idx, xs(idx));

            choice = questdlg(msg, ...
                'Manual polarity correction', ...
                'Flip and refit', 'Cancel', 'Flip and refit');

            if ~strcmp(choice, 'Flip and refit')
                continue
            end

            z_old = z_used;
            p_old = p;
            ci_old = ci;
            gof_old = gof;
            xs_old = xs;
            residual_old = residual_complex;
            flips_old = manual_flip_idx;

            z_used(idx) = -z_used(idx);

            R1start_refit = p(3);
            if ~isfinite(R1start_refit) || R1start_refit <= 0
                R1start_refit = R1start_centre;
            end

            [p2, ci2, gof2, ok2, xs2, z2, residual2] = fitComplexFcModelAuto( ...
                xs .* 1000, z_used, R1start_refit, lb, ub, lsqOpts, p);

            if ok2
                p = p2;
                ci = ci2;
                gof = gof2;
                xs = xs2;
                z_used = z2;
                residual_complex = residual2;

                if any(manual_flip_idx == idx)
                    manual_flip_idx(manual_flip_idx == idx) = [];
                else
                    manual_flip_idx(end+1) = idx;
                end
            else
                z_used = z_old;
                p = p_old;
                ci = ci_old;
                gof = gof_old;
                xs = xs_old;
                residual_complex = residual_old;
                manual_flip_idx = flips_old;

                warndlg('The polarity flip caused the refit to fail, so it was reverted.', ...
                    'Refit failed');
            end
        end

        if isgraphics(h)
            close(h);
        end

    end

    % ==========================================================
    % Nested helper: draw interactive fit figure
    % ==========================================================
    function drawManualComplexFitFigure(h, field_index, xs, z_used, p, gof, manual_flip_idx, auto_flip_idx)

        if ~isgraphics(h)
            return
        end

        figure(h);
        clf(h);

        A    = p(1);
        B    = p(2);
        R1   = p(3);
        phi0 = p(4);
        phi1 = p(5);

        phase_fit = phi0 + phi1 .* xs(:);
        z_rot = z_used .* exp(-1i .* phase_fit);

        y_real = real(z_rot);
        y_imag = imag(z_rot);

        xfit = linspace(min(xs), max(xs), 300).';
        yfit = A .* exp(-xfit .* R1) + B;

        plot(xs, y_real, 'o', ...
            'DisplayName', 'Real data after fitted phase correction');
        hold on;

        plot(xs, y_imag, 'x', ...
            'DisplayName', 'Imaginary component after fitted phase correction');

        plot(xfit, yfit, '-', ...
            'DisplayName', 'Real exponential fit');

        yline(0, '--', 'DisplayName', 'Zero');

        if nargin >= 8 && ~isempty(auto_flip_idx)
            auto_flip_idx = auto_flip_idx(auto_flip_idx >= 1 & auto_flip_idx <= numel(xs));
            if ~isempty(auto_flip_idx)
                plot(xs(auto_flip_idx), y_real(auto_flip_idx), 'd', ...
                    'MarkerSize', 10, ...
                    'LineWidth', 1.5, ...
                    'DisplayName', 'Automatically flipped point(s)');
            end
        end

        if ~isempty(manual_flip_idx)
            manual_flip_idx = manual_flip_idx(manual_flip_idx >= 1 & manual_flip_idx <= numel(xs));

            if ~isempty(manual_flip_idx)
                plot(xs(manual_flip_idx), y_real(manual_flip_idx), 's', ...
                    'MarkerSize', 10, ...
                    'LineWidth', 1.5, ...
                    'DisplayName', 'Manually flipped point(s)');
            end
        end

        hold off;
        grid on;

        xlabel('Evolution time (s)');
        ylabel('Signal (AU)');

        title(sprintf(['Field %d: R1 = %.4g s^{-1}, T1 = %.4g ms, ', ...
                       'A = %.4g, B = %.4g, phi0 = %.2f rad, phi1 = %.2f rad/s, R^2 = %.4f\n', ...
                       'Left-click point to flip polarity and refit. Enter/right-click accepts. u undoes.'], ...
            field_index, R1, 1000 ./ R1, A, B, phi0, phi1, gof.rsquare));

        legend('Location', 'best');
        drawnow;

    end

    % ==========================================================
    % Nested helper: static fit figure
    % ==========================================================
    function showStaticFitFigure(field_index, p, gof, xs, z_used, auto_flip_idx)

        h = figure;

        drawManualComplexFitFigure(h, field_index, xs, z_used, p, gof, [], auto_flip_idx);

        try
            waitforbuttonpress;
            close(h);
        catch
        end

    end

    % ==========================================================
    % Nested helper: nearest plotted real point
    % ==========================================================
    function idx = nearestRealFitPoint(xs, z_used, p, xclick, yclick)

        idx = [];

        if isempty(xs) || isempty(z_used) || numel(p) < 5
            return
        end

        phi0 = p(4);
        phi1 = p(5);

        phase_fit = phi0 + phi1 .* xs(:);
        z_rot = z_used .* exp(-1i .* phase_fit);

        y_real = real(z_rot);

        ax = gca;
        xl = xlim(ax);
        yl = ylim(ax);

        xrange = diff(xl);
        yrange = diff(yl);

        if xrange == 0 || ~isfinite(xrange)
            xrange = max(range(xs), eps);
        end

        if yrange == 0 || ~isfinite(yrange)
            yrange = max(range(y_real), eps);
        end

        d2 = ((xs(:) - xclick) ./ xrange).^2 + ...
             ((y_real(:) - yclick) ./ yrange).^2;

        [~, idx] = min(d2);

    end

    % ==========================================================
    % Nested helper: optional single-point trimming
    % ==========================================================
    function [bestP, bestCi, bestGof, okBest, bestXs, bestZ, bestResidual, bestRemoved] = ...
        fitComplexFcModelWithSingleTrim( ...
            x_ms, z, R1start_centre, lb, ub, lsqOpts, p_hint)

        okBest = false;

        bestP = nan(1, numel(lb));
        bestCi = nan(2, numel(lb));
        bestGof = emptyGof(Inf);
        bestXs = [];
        bestZ = [];
        bestResidual = [];
        bestRemoved = [];

        x_ms = x_ms(:);
        z = z(:);

        [x_ms_sorted, ord] = sort(x_ms);
        z_sorted = z(ord);

        nPts = numel(z_sorted);

        if nPts < 6
            return
        end

        for remove_idx = 1:nPts

            keep = true(nPts,1);
            keep(remove_idx) = false;

            [p_try, ci_try, gof_try, ok_try, xs_try, z_try, residual_try] = ...
                fitComplexFcModelAuto( ...
                    x_ms_sorted(keep), z_sorted(keep), ...
                    R1start_centre, lb, ub, lsqOpts, p_hint);

            if ~ok_try
                continue
            end

            if ~okBest || gof_try.rmse < bestGof.rmse
                bestP = p_try;
                bestCi = ci_try;
                bestGof = gof_try;
                bestXs = xs_try;
                bestZ = z_try;
                bestResidual = residual_try;
                bestRemoved = remove_idx;
                okBest = true;
            end
        end

    end

    % ==========================================================
    % Nested helper: approximate confidence intervals
    % ==========================================================
    function ci = approximateParameterCi(p, resnorm, J, dfe)

        ci = nan(2, numel(p));

        if isempty(J) || dfe <= 0 || ~all(isfinite(p))
            return
        end

        try
            mse = resnorm ./ dfe;
            JTJ = full(J.' * J);
            covp = mse .* pinv(JTJ);

            se = sqrt(max(diag(covp), 0)).';

            ci(1,:) = p - 1.96 .* se;
            ci(2,:) = p + 1.96 .* se;

        catch
            ci = nan(2, numel(p));
        end

    end

    % ==========================================================
    % Nested helper: empty gof struct
    % ==========================================================
    function gof = emptyGof(sseValue)

        gof = struct( ...
            'sse',        sseValue, ...
            'rsquare',    NaN, ...
            'adjrsquare', NaN, ...
            'dfe',        NaN, ...
            'rmse',       NaN, ...
            'exitflag',   NaN);

    end

    % ==========================================================
    % Nested helper: normalise parameter vector to [A B R1 phi0 phi1]
    % ==========================================================
    function p = normaliseParamVector(p)

        if isempty(p)
            return
        end

        p = p(:).';

        if numel(p) == 4
            p = [p, 0];
        elseif numel(p) > 5
            p = p(1:5);
        elseif numel(p) < 4
            p = [];
        end

    end

    % ==========================================================
    % Nested helper: wrap angle to [-pi, pi]
    % ==========================================================
    function phi = wrapPiLocal(phi)
        phi = mod(phi + pi, 2*pi) - pi;
    end

    % ==========================================================
    % Nested helper: property/field detection
    % ==========================================================
    function tf = hasProp(s, name)

        if isobject(s)
            tf = isprop(s, name);
        else
            tf = isfield(s, name);
        end

    end

    % ==========================================================
    % Nested helper: safe object/struct assignment
    % ==========================================================
    function s = safeSet(s, name, value)

        try
            if isobject(s)
                if isprop(s, name)
                    s.(name) = value;
                end
            else
                s.(name) = value;
            end
        catch
            % Do not fail the fit because a convenience property does not exist.
        end

    end

end

    end


end


