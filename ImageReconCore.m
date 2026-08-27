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
        fov_phase %phase-encoding FOV in mm
        fov_3d %through-plane/slab FOV in mm for 3D acquisitions
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
        fft_size_3d
        denoise_filter  %denoising filter
        denoise_params  %filter kernel
        mask
        dispersioncurve
        T1FitResults %complete independent/joint dispersion fit diagnostics
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
        fov_offsets_applied %true only for the current working k-space copy
        fov_offset_info %offset/FOV diagnostics for the most recent reconstruction
        selected_coils_original
        noise_whitening %apply receive-noise prewhitening during reconstruction
        whitening_info %diagnostics returned by noise_whiten
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
                    if isfield(fid.par.cameleon, 'FIELD_OF_VIEW_PHASE') && ...
                            isfinite(fid.par.cameleon.FIELD_OF_VIEW_PHASE) && ...
                            fid.par.cameleon.FIELD_OF_VIEW_PHASE > 0
                        obj.fov_phase = fid.par.cameleon.FIELD_OF_VIEW_PHASE*1000;
                    else
                        obj.fov_phase = obj.fov;
                    end
                    if isfield(fid.par.cameleon, 'FIELD_OF_VIEW_3D') && ...
                            isfinite(fid.par.cameleon.FIELD_OF_VIEW_3D) && ...
                            fid.par.cameleon.FIELD_OF_VIEW_3D > 0
                        obj.fov_3d = fid.par.cameleon.FIELD_OF_VIEW_3D*1000;
                    end
                    obj.TwoDimensional = fid.par.cameleon.MULTI_PLANAR_EXCITATION;
                    obj.thk = fid.par.cameleon.SLICE_THICKNESS*1000;
                    if isempty(obj.fov_3d)
                        obj.fov_3d = obj.thk*max(double(obj.slices), 1);
                    end
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
                    obj.fft_size_3d = obj.slices;
                    obj.denoise_filter = 'none';
                    obj.denoise_params = []; %filter kernel
                    obj.param = fid.par.cameleon;

                    if fid.par.cameleon.MULTI_PLANAR_EXCITATION ==1
                        obj.TwoDimensional = 1;
                    else
                        obj.TwoDimensional = 0; %ie 3d
                        obj.resolution_throughplane = obj.fov_3d/max(double(obj.slices), 1);
                        if isfield(fid.par.cameleon, 'USER_MATRIX_DIMENSION_3D')
                            userPartitions = round(double( ...
                                fid.par.cameleon.USER_MATRIX_DIMENSION_3D));
                            if isfinite(userPartitions) && userPartitions >= obj.slices
                                obj.fft_size_3d = userPartitions;
                            end
                        end
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

            % Backwards-compatible defaults shared by every input format.
            if isempty(obj.fov_phase)
                obj.fov_phase = obj.fov;
            end
            if isempty(obj.fov_3d)
                obj.fov_3d = obj.thk*max(double(obj.slices), 1);
            end
            hasMultipleReceivers = ~isempty(obj.n_receivers) && ...
                isscalar(obj.n_receivers) && isfinite(double(obj.n_receivers)) && ...
                double(obj.n_receivers) > 1;
            if isempty(obj.noise_whitening)
                obj.noise_whitening = hasMultipleReceivers;
            else
                obj.noise_whitening = logical(obj.noise_whitening) && hasMultipleReceivers;
            end
            obj.whitening_info = struct('method', 'notRun', 'nSamples', 0);
            obj.fov_offsets_applied = false;
            obj.fov_offset_info = struct('applied', false, ...
                'offsets_m', [0 0 0], 'fov_m', [NaN NaN NaN], ...
                'is3D', false);
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
                % Keep an untouched, unwhitened copy. Whitening is a user-selectable
                % reconstruction operation and must be applied exactly once to
                % the full receiver array in buildimages(), before selection.
                obj.originalcomplexkspace = obj.complexkspace;
               
            else
                obj.complexkspace = reshape( obj.complexkspace,[obj.samples,size(obj.complexkspace,2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]); %enforce dimensionality
                obj.originalcomplexkspace = obj.complexkspace;
            end

            % originalcomplexkspace is always the immutable, unshifted source.
            % FOV offsets are applied only to a fresh working copy in
            % buildimages().
            obj.fov_offsets_applied = false;

        end


        function obj = buildimages(obj)

            %BUILDIMAGES  Convert k-space to image space, perform coil combination, optional bias correction, filtering.
            %
            
            obj.compleximage = [];
            % Always restart from immutable preprocessed acquisition data.
            % This makes FOV-offset correction independent of the number of
            % rebuilds and of any accumulated display rotation/flip.
            obj.complexkspace = obj.originalcomplexkspace;
            obj.fov_offsets_applied = false;
            obj = correct_orientation(obj);

            % ----------- Receiver selection and optional prewhitening -----------
            % obj.multichannel_recon may be:
            %   - logical mask over the ORIGINAL coil dimension
            %   - numeric MATLAB index list into the ORIGINAL coil dimension
            %     (the GUI displays scanner receiver IDs 0..N-1, with
            %     ItemsData 1..N)
            opts_orig = obj.multichannel_recon;
            nCoils_orig = size(obj.complexkspace, 6);

            if isempty(opts_orig)
                sel_orig = 1:nCoils_orig;
            elseif islogical(opts_orig)
                if numel(opts_orig) ~= nCoils_orig
                    error('buildimages:coilMaskWrongLength', ...
                        'Logical receiver mask has %d entries; %d receivers are available.', ...
                        numel(opts_orig), nCoils_orig);
                end
                sel_orig = find(opts_orig(:).');
            else
                if ~isnumeric(opts_orig)
                    error('buildimages:coilSelNotNumeric', ...
                        'multichannel_recon must be a numeric index list or logical mask. Got %s.', ...
                        class(opts_orig));
                end

                sel_orig = unique(round(opts_orig(:).'), 'stable');
                sel_orig(~isfinite(sel_orig)) = [];

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
            nSel = numel(sel_orig);

            % Estimate and apply one whitening transform across the complete
            % acquired receiver array.  Selection is deliberately applied
            % afterwards: for example, selecting displayed receiver 7 retains
            % MATLAB channel 8 of the fully prewhitened array.  A retained
            % channel is therefore a virtual (mixed) prewhitened channel, not
            % the untouched physical receiver with that number.
            whiteningEnabled = ~isempty(obj.noise_whitening) && ...
                isscalar(obj.noise_whitening) && logical(obj.noise_whitening) && ...
                nCoils_orig > 1;
            [kspace_selected, obj.whitening_info] = prewhiten_and_select( ...
                obj.complexkspace, sel_orig, whiteningEnabled);
            opts = 1:nSel;

            % Magnitude k-space is kept at the acquired matrix size for QC.
            obj.magkspace = rssq(kspace_selected, 6);

            % ----------- Precompute centred zero filling -----------
            % Spatial dimension order is read / phase / partition.
            nRead  = size(kspace_selected, 1);
            nPhase = size(kspace_selected, 2);
            nPart  = size(kspace_selected, 3);

            targetRead = round(double(obj.fft_size));
            if ~isfinite(targetRead) || targetRead < 1
                error('buildimages:invalidFftSize', ...
                    'obj.fft_size must be a positive finite integer.');
            end

            % Apply the same zero-fill factor to read and phase. Deriving this
            % from the actual array sizes avoids stale header values and works
            % for rectangular matrices/non-square FOV acquisitions.
            targetPhase = round(double(nPhase) * targetRead / double(nRead));
            targetPart = nPart;
            if obj.TwoDimensional ~= 1 && ~isempty(obj.fft_size_3d)
                targetPart = round(double(obj.fft_size_3d));
            end

            if targetRead < nRead
                error('buildimages:fftSizeTooSmall', ...
                    'obj.fft_size=%d is smaller than acquired readout size %d.', ...
                    targetRead, nRead);
            end
            if targetPhase < nPhase
                error('buildimages:phaseFftSizeTooSmall', ...
                    'Target phase FFT size %d is smaller than acquired phase size %d.', ...
                    targetPhase, nPhase);
            end
            if ~isfinite(targetPart) || targetPart < nPart
                error('buildimages:partitionFftSizeTooSmall', ...
                    'obj.fft_size_3d=%g is smaller than acquired partition size %d.', ...
                    double(targetPart), nPart);
            end

            padReadTotal  = targetRead  - nRead;
            padPhaseTotal = targetPhase - nPhase;
            padPartTotal  = targetPart  - nPart;

            % padarray entries follow the real array order: read, phase, part.
            % The previous phase/read reversal only went unnoticed for square data.
            padPre = [floor(padReadTotal/2), ...
                floor(padPhaseTotal/2), floor(padPartTotal/2)];
            padPost = [ceil(padReadTotal/2), ...
                ceil(padPhaseTotal/2), ceil(padPartTotal/2)];

            kspace_selected = windowkspace(kspace_selected, ...
                obj.window_size, obj.window_function);
            kpad = padarray(kspace_selected, padPre, 0, 'pre');
            kpad = padarray(kpad, padPost, 0, 'post');

            % ----------- FFT to image space (per coil) -----------
            if obj.TwoDimensional == 1
                obj.compleximage = ifft2c(kpad);
            else
                % Transform only the first three spatial dimensions. The old
                % ifftn implementation also transformed the receiver dimension.
                obj.compleximage = ifft3c(kpad);
            end

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
                        ccopts = struct();
                        ccopts.coils = 1:size(obj.compleximage, 6);
                        ccopts.sigma = 4;
                        ccopts.maskFrac = 0.1;
                        if whiteningEnabled
                            % Whitening creates virtual channels; no output
                            % channel remains the privileged "signal coil".
                            ccopts.referenceMode = 'rss';
                        else
                            ccopts.referenceMode = 'signalcoil';
                        end
                        obj.complexcombined = combine_adaptive_smoothed(obj.compleximage, ccopts);
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
                        % Bias correction previously referenced `maps` even
                        % when RSSQ or Adaptive combination had not created it.
                        if ~exist('maps', 'var') || isempty(maps)
                            maps = estimate_maps_from_acs_walsh(kpad, opts, obj);
                        end
                        obj.magimage = map_guided_bias_correct(obj.magimage, maps, obj);

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

            % Apply the accumulated lossless display orientation only after
            % offset correction, FFT and receiver combination. Raw k-space is
            % never rotated or interpolated.
            if ~isempty(obj.geomT)
                obj = apply_reconstruction_geometry(obj, obj.geomT, false);
            end
        end


     function obj = MapRelaxation(obj)
%MAPRELAXATION Pixelwise FCI T1 mapping, parfor-safe class wrapper.
%
% This wrapper intentionally contains NO parfor loop.  The parallel work is
% done in the standalone function ffc_t1map_phase_signed_run_parallel.m so
% MATLAB workers do not need to execute ImageReconCore.m or any local helper
% functions inside the class file.
%
% Required files on the MATLAB path:
%   ffc_t1map_phase_signed_run_parallel.m
%   ffc_t1map_phase_signed_magnitude_worker_v2.m
%
% Model per pixel/field/slice:
%   y(t) = B - C*exp(-R1*t), B>=0, C>=0, R1>=0
%
% If compleximage is available, phase is used only to assign sign while
% preserving magnitude:
%   y = sign(real(z*exp(-1i*phi_ref))) .* abs(z)

    opts = struct();
    opts.MaskFactor = 0.02;
    opts.UseExistingMask = false;
    opts.PhaseReferenceMode = 'last2';
    opts.SignAssignmentMode = 'phase_then_prefix_if_bad';
    opts.MaxPrefixSignCandidates = 4;
    opts.T1AcceptedRangeMs = [50 3000];
    opts.RsqBadThreshold = 0.85;
    opts.R1Min = 0.2;
    opts.R1Max = 40;
    opts.R1GridPoints = 80;
    opts.R1GridLogSpaced = true;
    opts.MinValidTimepoints = 4;
    opts.UseParallel = true;

    % Optional overrides can be placed in obj.param.T1MapOptions as a struct.
    try
        if isstruct(obj.param) && isfield(obj.param, 'T1MapOptions') && isstruct(obj.param.T1MapOptions)
            userOpts = obj.param.T1MapOptions;
            fn = fieldnames(userOpts);
            for k = 1:numel(fn)
                opts.(fn{k}) = userOpts.(fn{k});
            end
        end
    catch
    end

    if isempty(obj.magimage)
        error('MapRelaxation:NoMagnitudeImage', 'obj.magimage is empty. Rebuild images before T1 mapping.');
    end

    magData = obj.magimage;

    useComplex = false;
    complexData = [];
    try
        if ~isempty(obj.compleximage) && ndims(obj.compleximage) >= 4 && ...
                size(obj.compleximage,1) == size(magData,1) && size(obj.compleximage,2) == size(magData,2)
            complexData = obj.compleximage;
            useComplex = true;
        end
    catch
        useComplex = false;
        complexData = [];
    end

    times = double(obj.timepoints);
    fields = double(obj.fieldpoints(:)); %#ok<NASGU>

    if isscalar(obj.n_fieldpoints)
        nFields = double(obj.n_fieldpoints);
    else
        nFields = double(obj.n_fieldpoints(1));
    end

    % A 3D reconstruction may be zero-filled through-plane, so the output
    % stack can contain more planes than the acquired partition count.
    nSlices = size(magData, 3);

    userMask = [];
    try
        if ~isempty(obj.mask)
            userMask = logical(obj.mask);
        end
    catch
        userMask = [];
    end

    out = ffc_t1map_phase_signed_run_parallel( ...
        magData, complexData, useComplex, times, nFields, nSlices, userMask, opts);

    obj.T1Maps = out.T1;
    obj.R1Maps = out.R1;

    % Store extra QC outputs in obj.param so you do not have to add new class
    % properties.  Add real properties later if you want direct dot access.
    try
        obj.param.T1Map_Rsq = out.Rsq;
        obj.param.T1Map_RMSE = out.RMSE;
        obj.param.T1Map_SSE = out.SSE;
        obj.param.T1Map_QC = out.QC;
        obj.param.T1Map_SignMode = out.SignMode;
        obj.param.T1Map_Mask = out.Mask;
        obj.param.T1MapOptionsUsed = opts;
    catch
    end

    fprintf('Pixelwise T1 mapping complete: %d x %d x %d slices x %d fields. Complex phase-signing: %d\n', ...
        size(out.T1,1), size(out.T1,2), size(out.T1,3), size(out.T1,4), useComplex);
end



function [obj, results] = T1dispersion(obj, slice, userOptions)
%T1DISPERSION Phase-aware field-cycling inversion-recovery T1 fit.
%
%   OBJ = T1DISPERSION(OBJ, SLICE)
%   [OBJ, RESULTS] = T1DISPERSION(OBJ, SLICE, OPTIONS)
%
% This is a drop-in replacement for the previous field-wise fitter.  It is
% designed for inversion-recovery field-cycling experiments in which the
% equilibrium magnetisation changes during the evolution field and again
% during the return to the detection field.
%
% If all preparation, ramp and readout timings other than the evolution time
% are fixed, the longitudinal signal at one evolution field can always be
% written as
%
%       S(t) = S_inf + (S_0 - S_inf) exp(-R1 t)
%            = B - C exp(-R1 t).                         (1)
%
% B and C are nuisance parameters and are NOT tied to the detection-field
% M0. In the default inversion mode B>=0 and C>B, so the initial state is
% negative while a near-zero equilibrium remains possible at ultra-low field.
% Set sequenceMode='unconstrained' only for non-IR acquisitions.
%
% Primary fitting routes
% ---------------------
% The default auto mode first tests whether a complex ROI signal can be
% described by a smoothly drifting global phase plus a single IR sign
% transition.  If so, the drift-corrected signed complex signal is fitted.
% Otherwise, the magnitude is fitted directly without assigning measured
% samples arbitrary signs.  For a single T/R coil the expected Rician
% magnitude is used; for combined multi-coil magnitude data a measured-noise
% floor approximation is used.  Both routes enforce an inverted initial
% state, S(0)<0, unless OPTIONS.sequenceMode is explicitly changed.
%
% Optional joint multi-field fit
% ------------------------------
% Set OPTIONS.fitMode='joint' to refine phase-validated signed field-wise
% fits using the combined-field signal model of Boedenler et al.:
%
%   S_f(t) = M0_D*[r_f*(1-exp(-R1_f*t))
%                         - alpha_f*exp(-R1_f*t)],          (2)
%
% where r_f=B_E,f/B_D, M0_D is shared by every field, and alpha_f and R1_f
% are field-specific.  In an inversion-recovery sequence alpha_f must be
% non-negative: S_f(0)=-M0_D*alpha_f.  The linear subproblem therefore uses
% A_f=M0_D*(alpha_f-alpha_min)>=0 rather than treating the exponential
% amplitude as an unconstrained nuisance parameter.  This prevents a
% magnitude curve from being explained by an unphysical negative inversion
% efficiency.  The shared equilibrium scaling is what makes this a genuinely
% simultaneous fit; merely summing independent field objectives would return
% exactly the original field-wise estimates.
%
% The independent mode is the default because the joint model additionally
% requires intensity calibration consistent with M0 proportional to field.
% This ROI implementation reproduces the paper's combined-field pixel-wise model;
% it does not implement the paper's k-space Frobenius-TGV reconstruction.
% Reference: Boedenler M et al., Magn Reson Med. 2021;86:2049-2063,
% doi:10.1002/mrm.28857.
%
% Automatic complex/magnitude fit
% -------------------------------
% In the default 'auto' mode, every physically possible IR zero-crossing
% branch is tested.  For each branch a low-order, smoothly varying phase drift
% is estimated without independently rotating the time points.  This preserves
% a genuine pi inversion while allowing ppm-level readout-field drift.
%
% Magnitude inversion-recovery model
% ----------------------------------
% Magnitude observations are fitted as a function of
%       nu(t)=abs(B-C*exp(-R1*t)),  B>=0, C>B,
% rather than converted into signed Gaussian samples.  This explicitly models
% the non-zero magnitude near the inversion null.  OPTIONS.magnitudeNoiseSigma
% can provide a measured single-quadrature noise SD; otherwise a background
% estimate is used and the recovery minimum is only a last-resort estimate.
%
% Robustness and uncertainty
% --------------------------
% Huber iteratively reweighted least squares reduces the influence of motion
% or corrupted time points.  A one-dimensional profile-likelihood interval
% is calculated for R1; this is more informative than a symmetric Jacobian
% interval when the asymptote is close to zero or R1 approaches a bound.
%
% Required OBJ members
% --------------------
%   fieldpoints    evolution fields in mT
%   timepoints     evolution times in ms; vector or [field x time]
%   mask           2-D mask or [X Y slice] mask
%   compleximage   preferred, [X Y slice time field]
%       OR
%   magimage       fallback, [X Y slice time field]
%
% Existing output members are populated when they exist on an object, and
% are always added when OBJ is a struct.  Detailed diagnostics are returned
% in RESULTS even when the OBJ class has no matching property.
%
% Useful OPTIONS fields (all optional)
% ------------------------------------
%   signalMode          'auto' (default), 'complex', or 'magnitude'.
%                       'auto' fits both representations and uses complex
%                       results only when their phase/polarity is coherent.
%   roiAggregation      'phase_aligned_trimmed_mean' (default),
%                       'phase_aligned_mean', or 'plain_mean'
%   trimFraction        spatial trimming fraction, default 0.10
%   robust              false (default); with few time points robust and
%                       ordinary fits should be compared rather than silently
%                       discarding influential recovery samples
%   robustTune          radial Huber tuning constant, default 2.5
%   robustIterations    default 8
%   T1BoundsMs          [20 3000] by default
%   fitMode             'independent' (default) or 'joint'
%   sequenceMode        'inversion' (default) or 'unconstrained'
%   minimumInitialMagnetizationFraction  minimum |S(0)| relative to the
%                       observed signal scale, default 0.05
%   phaseDriftOrder     polynomial order versus stored acquisition index,
%                       default 1
%   maximumPhaseResidualDegrees  maximum RMS residual for accepting complex
%                       phase, default 40 degrees
%   magnitudeNoiseSigma NaN (auto), scalar, or one value per field
%   magnitudeCoilCount NaN (infer), or number of coils represented by a
%                       magnitude-only image
%   forcedNegativeIndices / forcedPositiveIndices  cell arrays indexed by
%                       field, used to impose known sample polarity after
%                       evolution-time sorting and duplicate-time averaging
%   detectionField_mT   NaN (default): infer max(fieldpoints). Set the true
%                       detection/inversion field explicitly when it is not
%                       the largest evolution field; it enters the physical
%                       inversion constraint as well as reported alpha.
%   jointRateSweeps     bounded coordinate sweeps per robust iteration,
%                       default 5
%   jointSignIterations joint magnitude-polarity refinements, default 3
%   jointMinimumAlpha   minimum inversion/ramp efficiency, default 0.10.
%                       Set to 0 only if a non-inverted initial state is
%                       physically possible in the acquisition.
%   jointTimeWeighting  false (default): optional later-time upweighting
%                       evolution times within each field
%   jointLateTimeWeight relative latest/earliest weight, default 4
%   jointTimeWeightExponent curvature of the normalized time ramp,
%                       default 1 (linear)
%   jointNoiseWeighting weight fields by robust residual noise, default true
%   makeDispersionPlot  true (default)
%   showDashboard       defaults to obj.checkfit == 1
%   printSummary        true (default)
%
% R1error and T1error are approximate one-standard-error values derived from
% the profile interval.  The full intervals are in results.field(n).R1CI95
% and results.field(n).T1CI95, and in obj.R1ci95 / obj.T1ci95 when supported.

    if nargin < 2 || isempty(slice)
        slice = 1;
    end
    if nargin < 3 || isempty(userOptions)
        userOptions = struct();
    end
    if ~isstruct(userOptions)
        error('T1dispersion:InvalidOptions', 'OPTIONS must be a struct.');
    end

    opts = defaultOptions();
    opts = mergeOptions(opts, userOptions);
    if ~isfield(userOptions, 'showDashboard')
        opts.showDashboard = memberExists(obj, 'checkfit') && ...
            isscalar(obj.checkfit) && obj.checkfit == 1;
    end
    opts = validateOptions(opts);

    if ~memberExists(obj, 'fieldpoints') || isempty(obj.fieldpoints)
        error('T1dispersion:MissingFields', 'obj.fieldpoints is required.');
    end
    if ~memberExists(obj, 'timepoints') || isempty(obj.timepoints)
        error('T1dispersion:MissingTimes', 'obj.timepoints is required.');
    end
    if ~memberExists(obj, 'mask') || isempty(obj.mask)
        error('T1dispersion:MissingMask', 'obj.mask is required.');
    end

    fields_mT = double(obj.fieldpoints(:));
    times_ms = double(obj.timepoints);
    mask2d = selectMask(obj.mask, slice);
    if ~any(mask2d(:))
        error('T1dispersion:EmptyMask', 'The selected mask is empty for slice %d.', slice);
    end

    % fieldpoints is the authoritative list.  Do not use the currently
    % displayed EvolutionFieldDropDown value and do not silently reduce the
    % fit to a four-dimensional, displayed-field image.
    nFields = numel(fields_mT);
    sources = resolveSignalSources(obj, opts.signalMode, nFields, mask2d, slice);
    if isnan(opts.magnitudeCoilCount)
        opts.magnitudeCoilCount = inferMagnitudeCoilCount(obj, sources);
    end
    signalMode = sources.mode;
    implementationTag = ['ImageReconCore.T1dispersion parent-integrated ', ...
        'v8 phase-drift-rician-constrained-IR'];
    inputDiagnostics = struct( ...
        'fieldCount', nFields, ...
        'timepointsSize', size(times_ms), ...
        'magimageSize', memberSize(obj, 'magimage'), ...
        'compleximageSize', memberSize(obj, 'compleximage'), ...
        'selectedMode', signalMode, ...
        'fitMode', opts.fitMode, ...
        'complexMember', sources.complexMember, ...
        'magnitudeMember', sources.magnitudeMember, ...
        'complexRejectedReason', sources.complexRejectedReason);
    inputDiagnostics.magnitudeCoilCount = opts.magnitudeCoilCount;

    if memberExists(obj, 'n_fieldpoints') && ~isempty(obj.n_fieldpoints)
        declaredCount = double(obj.n_fieldpoints(1));
        if isfinite(declaredCount) && declaredCount ~= nFields
            warning('T1dispersion:DeclaredFieldCountMismatch', ...
                ['obj.n_fieldpoints is %d but obj.fieldpoints contains %d ', ...
                 'values. Fitting all %d values in obj.fieldpoints.'], ...
                declaredCount, nFields, nFields);
        end
    end

    fieldTemplate = emptyFieldResult();
    fieldResult = repmat(fieldTemplate, 1, nFields);

    % This method contains nested helper functions. MATLAB shares a parent
    % variable with a nested helper when both use the same name. The helpers
    % use `n` for the number of time samples, so the field loop must use a
    % unique name or every fit is written into the same result slot.
    for fieldLoopIndex = 1:nFields
        try
            t_ms = timeVectorForField(times_ms, fieldLoopIndex, nFields);
            complexSignal = [];
            magnitudeSignal = [];
            complexMeta = struct();
            magnitudeMeta = struct();

            if sources.useComplex
                [complexSignal, complexMeta] = extractRoiSignal( ...
                    sources.complexData, mask2d, slice, fieldLoopIndex, ...
                    'complex', opts);
            end
            if sources.useMagnitude
                [magnitudeSignal, magnitudeMeta] = extractRoiSignal( ...
                    sources.magnitudeData, mask2d, slice, fieldLoopIndex, ...
                    'magnitude', opts);
            elseif ~isempty(complexSignal)
                % Magnitude is taken only after complex ROI averaging, which
                % avoids per-voxel Rician bias as far as the available data allow.
                magnitudeSignal = abs(complexSignal);
                magnitudeMeta = complexMeta;
                magnitudeMeta.aggregation = 'magnitude_of_complex_roi_mean';
                if isfield(magnitudeMeta, 'roiNoiseSigma') && ...
                        isfinite(magnitudeMeta.roiNoiseSigma)
                    magnitudeMeta.backgroundMagnitudeMean = ...
                        magnitudeMeta.roiNoiseSigma .* sqrt(pi ./ 2);
                end
            end

            [complexTime, complexSignal] = prepareFieldSamples( ...
                t_ms, complexSignal, fieldLoopIndex, 'complex');
            [magnitudeTime, magnitudeSignal] = prepareFieldSamples( ...
                t_ms, magnitudeSignal, fieldLoopIndex, 'magnitude');

            switch signalMode
                case 'complex'
                    polarity = fieldPolarityConstraints(obj, opts, ...
                        fieldLoopIndex, numel(complexTime));
                    currentFieldResult = fitOrFailComplex(complexTime, complexSignal, ...
                        opts, fieldTemplate, polarity, fieldLoopIndex);
                    currentFieldResult.roiMeta = complexMeta;
                case 'magnitude'
                    polarity = fieldPolarityConstraints(obj, opts, ...
                        fieldLoopIndex, numel(magnitudeTime));
                    currentFieldResult = fitOrFailMagnitude( ...
                        magnitudeTime, magnitudeSignal, ...
                        opts, fieldTemplate, polarity, fieldLoopIndex, ...
                        magnitudeMeta);
                    currentFieldResult.roiMeta = magnitudeMeta;
                otherwise % hybrid_auto
                    polarity = fieldPolarityConstraints(obj, opts, ...
                        fieldLoopIndex, max(numel(complexTime), ...
                        numel(magnitudeTime)));
                    complexFit = fitOrFailComplex(complexTime, complexSignal, ...
                        opts, fieldTemplate, polarity, fieldLoopIndex);
                    magnitudeFit = fitOrFailMagnitude(magnitudeTime, magnitudeSignal, ...
                        opts, fieldTemplate, polarity, fieldLoopIndex, ...
                        magnitudeMeta);
                    currentFieldResult = selectHybridFit( ...
                        complexFit, magnitudeFit, opts);
                    currentFieldResult.roiMeta = struct('complex', complexMeta, ...
                        'magnitude', magnitudeMeta);
            end
            currentFieldResult.fieldIndex = fieldLoopIndex;
            currentFieldResult.field_mT = fields_mT(fieldLoopIndex);
            fieldResult(fieldLoopIndex) = currentFieldResult;
        catch fieldException
            currentFieldResult = fieldTemplate;
            currentFieldResult.fieldIndex = fieldLoopIndex;
            currentFieldResult.field_mT = fields_mT(fieldLoopIndex);
            currentFieldResult.method = signalMode;
            currentFieldResult.qc = ['failed: ' fieldException.message];
            currentFieldResult.qcReasons = ...
                {fieldException.identifier, fieldException.message};
            fieldResult(fieldLoopIndex) = currentFieldResult;
            warning('T1dispersion:FieldFitFailed', ...
                'Field %d (%.5g mT) failed: %s', fieldLoopIndex, ...
                fields_mT(fieldLoopIndex), fieldException.message);
        end
    end

    % Every field slot must have been touched, whether the individual fit
    % succeeded or failed. This turns traversal corruption into an immediate
    % diagnostic instead of silently printing nine "not fitted" rows.
    untouchedFieldSlots = find(strcmp({fieldResult.qc}, 'not fitted'));
    writtenFieldIndices = [fieldResult.fieldIndex];
    if ~isempty(untouchedFieldSlots) || ...
            ~isequal(writtenFieldIndices, 1:nFields)
        error('T1dispersion:InternalFieldTraversalFailure', ...
            ['The field loop did not populate every result slot. Untouched ', ...
             'slots: %s; written indices: %s.'], ...
            mat2str(untouchedFieldSlots), mat2str(writtenFieldIndices));
    end

    independentFieldResult = fieldResult;
    jointDiagnostics = emptyJointDiagnostics(fields_mT, opts);
    if strcmp(opts.fitMode, 'joint')
        try
            [fieldResult, jointDiagnostics] = fitJointMultiField( ...
                fieldResult, fields_mT, opts);
        catch jointException
            jointDiagnostics.attempted = true;
            jointDiagnostics.used = false;
            jointDiagnostics.failureIdentifier = jointException.identifier;
            jointDiagnostics.failureMessage = jointException.message;
            if opts.jointFallbackToIndependent
                warning('T1dispersion:JointFitFailed', ...
                    ['Joint multi-field fitting failed; returning the ', ...
                     'independent initial fits. Cause: %s'], ...
                    jointException.message);
                fieldResult = independentFieldResult;
            else
                rethrow(jointException);
            end
        end
    end

    results = assembleResults(fieldResult, fields_mT, slice, signalMode, opts);
    results.implementation = implementationTag;
    results.inputDiagnostics = inputDiagnostics;
    results.joint = jointDiagnostics;
    results.independentInitial = arrayfun(@fitSummary, ...
        independentFieldResult, 'UniformOutput', false);
    if jointDiagnostics.used
        results.modelEquation = jointDiagnostics.modelEquation;
        results.modelName = ...
            'joint_multifield_inversion_constrained_variable_projection';
    end
    obj = storeResults(obj, results);

    if opts.printSummary
        printFitSummary(results);
    end
    if opts.makeDispersionPlot
        makeDispersionPlot(results, obj);
    end
    if opts.showDashboard
        makeQaDashboard(results);
    end


% ========================================================================
% Nested helpers
%
% These functions deliberately remain inside T1dispersion rather than being
% file-local subfunctions. ImageReconCore defines T1dispersion as a class
% method; in that context file-local helpers may instead be parsed as class
% methods and a bare call such as defaultOptions() is then unresolved.
% Nesting keeps the implementation self-contained in both installation
% layouts.
% ========================================================================

% ========================================================================
% Options and input handling
% ========================================================================

function opts = defaultOptions()
    opts = struct();
    opts.fitMode = 'independent';
    opts.signalMode = 'magnitude';
    opts.sequenceMode = 'inversion';
    opts.roiAggregation = 'phase_aligned_trimmed_mean';
    opts.trimFraction = 0.10;
    opts.anchorMinFraction = 0.05;
    opts.robust = false;
    opts.robustTune = 2.5;
    opts.robustIterations = 8;
    opts.robustTolerance = 1e-3;
    opts.minimumWeight = 0.05;
    opts.outlierWeightThreshold = 0.70;
    opts.T1BoundsMs = [20 3000];
    opts.rateGridPoints = 120;
    opts.profileGridPoints = 320;
    opts.confidenceLevel = 0.95;
    opts.minimumPoints = 6;
    opts.minimumInitialMagnetizationFraction = 0.05;
    opts.phaseDriftOrder = 1;
    opts.maximumPhaseResidualDegrees = 40;
    opts.minimumPhaseSNR = 3;
    opts.phaseReliableMagnitudeFraction = 0.15;
    opts.magnitudeNoiseSigma = NaN;
    opts.magnitudeNoiseModel = 'auto';
    opts.magnitudeCoilCount = NaN;
    opts.backgroundFraction = 0.20;
    opts.forcedNegativeIndices = {};
    opts.forcedPositiveIndices = {};
    opts.signAmbiguityRatio = 1.10;
    opts.signAmbiguityLogRateDifference = 0.15;
    opts.boundCandidatePenaltyFraction = 0.05;
    opts.maximumReliableSignTransitions = 1;
    opts.maximumSignIndexDisagreement = 1;
    opts.minimumDynamicSNR = 4;
    opts.maximumNRMSE = 0.20;
    opts.maximumOrthogonalNRMSE = 0.20;
    opts.minimumLateCoverageT1 = 1.5;
    opts.maximumEarlyCoverageT1 = 0.5;
    opts.maximumT1CIRatio = 2.0;
    opts.detectionField_mT = 200;
    opts.jointMinimumFields = 2;
    opts.jointRateSweeps = 5;
    opts.jointRobustIterations = 5;
    opts.jointSignIterations = 3;
    opts.jointTolerance = 1e-3;
    opts.jointMinimumAlpha = 0.10;
    opts.jointTimeWeighting = false;
    opts.jointLateTimeWeight = 4;
    opts.jointTimeWeightExponent = 1;
    opts.jointNoiseWeighting = true;
    opts.jointMaximumNoiseWeightRatio = 25;
    opts.jointFallbackToIndependent = true;
    opts.makeDispersionPlot = true;
    opts.showDashboard = true;
    opts.printSummary = true;
end

function opts = mergeOptions(opts, userOptions)
    names = fieldnames(userOptions);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(opts, name)
            error('T1dispersion:UnknownOption', 'Unknown option "%s".', name);
        end
        opts.(name) = userOptions.(name);
    end
end

function opts = validateOptions(opts)
    validFitModes = {'joint', 'independent'};
    if ~ischar(opts.fitMode) && ~(isstring(opts.fitMode) && isscalar(opts.fitMode))
        error('T1dispersion:InvalidFitMode', 'fitMode must be text.');
    end
    opts.fitMode = lower(char(opts.fitMode));
    if ~any(strcmp(opts.fitMode, validFitModes))
        error('T1dispersion:InvalidFitMode', ...
            'fitMode must be joint or independent.');
    end

    validSignalModes = {'auto', 'complex', 'magnitude'};
    if ~ischar(opts.signalMode) && ~(isstring(opts.signalMode) && isscalar(opts.signalMode))
        error('T1dispersion:InvalidSignalMode', 'signalMode must be text.');
    end
    opts.signalMode = lower(char(opts.signalMode));
    if ~any(strcmp(opts.signalMode, validSignalModes))
        error('T1dispersion:InvalidSignalMode', ...
            'signalMode must be auto, complex, or magnitude.');
    end

    validSequenceModes = {'inversion', 'unconstrained'};
    if ~ischar(opts.sequenceMode) && ...
            ~(isstring(opts.sequenceMode) && isscalar(opts.sequenceMode))
        error('T1dispersion:InvalidSequenceMode', ...
            'sequenceMode must be text.');
    end
    opts.sequenceMode = lower(char(opts.sequenceMode));
    if ~any(strcmp(opts.sequenceMode, validSequenceModes))
        error('T1dispersion:InvalidSequenceMode', ...
            'sequenceMode must be inversion or unconstrained.');
    end

    validMagnitudeModels = {'auto', 'rician', 'noncentral_chi_approx'};
    opts.magnitudeNoiseModel = lower(char(opts.magnitudeNoiseModel));
    if ~any(strcmp(opts.magnitudeNoiseModel, validMagnitudeModels))
        error('T1dispersion:InvalidMagnitudeNoiseModel', ...
            ['magnitudeNoiseModel must be auto, rician, or ', ...
             'noncentral_chi_approx.']);
    end

    validRoiModes = {'phase_aligned_trimmed_mean', 'phase_aligned_mean', 'plain_mean'};
    opts.roiAggregation = lower(char(opts.roiAggregation));
    if ~any(strcmp(opts.roiAggregation, validRoiModes))
        error('T1dispersion:InvalidRoiMode', ...
            'Unsupported roiAggregation "%s".', opts.roiAggregation);
    end
    if ~isscalar(opts.trimFraction) || opts.trimFraction < 0 || opts.trimFraction >= 0.5
        error('T1dispersion:InvalidTrim', 'trimFraction must be in [0, 0.5).');
    end
    if numel(opts.T1BoundsMs) ~= 2 || any(~isfinite(opts.T1BoundsMs)) || ...
            any(opts.T1BoundsMs <= 0) || opts.T1BoundsMs(1) == opts.T1BoundsMs(2)
        error('T1dispersion:InvalidBounds', ...
            'T1BoundsMs must contain two distinct positive finite values.');
    end
    opts.T1BoundsMs = sort(double(opts.T1BoundsMs(:).'));
    opts.R1Bounds = [1000 / opts.T1BoundsMs(2), 1000 / opts.T1BoundsMs(1)];

    scalarNonnegative = {'minimumInitialMagnetizationFraction', ...
        'maximumPhaseResidualDegrees'};
    for optionIndex = 1:numel(scalarNonnegative)
        optionName = scalarNonnegative{optionIndex};
        optionValue = opts.(optionName);
        if ~isscalar(optionValue) || ~isfinite(optionValue) || optionValue < 0
            error('T1dispersion:InvalidOption', ...
                '%s must be a non-negative finite scalar.', optionName);
        end
    end
    if ~isscalar(opts.phaseDriftOrder) || ~isfinite(opts.phaseDriftOrder) || ...
            opts.phaseDriftOrder < 0 || fix(opts.phaseDriftOrder) ~= opts.phaseDriftOrder
        error('T1dispersion:InvalidPhaseDriftOrder', ...
            'phaseDriftOrder must be a non-negative integer.');
    end
    if ~isscalar(opts.minimumPhaseSNR) || ~isfinite(opts.minimumPhaseSNR) || ...
            opts.minimumPhaseSNR <= 0
        error('T1dispersion:InvalidMinimumPhaseSNR', ...
            'minimumPhaseSNR must be a positive finite scalar.');
    end
    if ~isscalar(opts.phaseReliableMagnitudeFraction) || ...
            ~isfinite(opts.phaseReliableMagnitudeFraction) || ...
            opts.phaseReliableMagnitudeFraction < 0 || ...
            opts.phaseReliableMagnitudeFraction >= 1
        error('T1dispersion:InvalidPhaseMagnitudeFraction', ...
            'phaseReliableMagnitudeFraction must be in [0,1).');
    end
    if ~isscalar(opts.backgroundFraction) || ~isfinite(opts.backgroundFraction) || ...
            opts.backgroundFraction <= 0 || opts.backgroundFraction > 1
        error('T1dispersion:InvalidBackgroundFraction', ...
            'backgroundFraction must be in (0,1].');
    end
    if ~(isnumeric(opts.magnitudeNoiseSigma) && ...
            all(isfinite(opts.magnitudeNoiseSigma(:)) | ...
                isnan(opts.magnitudeNoiseSigma(:))) && ...
            all(opts.magnitudeNoiseSigma(isfinite(opts.magnitudeNoiseSigma)) > 0))
        error('T1dispersion:InvalidMagnitudeNoiseSigma', ...
            'magnitudeNoiseSigma must contain NaN or positive values.');
    end
    if ~(isscalar(opts.magnitudeCoilCount) && ...
            (isnan(opts.magnitudeCoilCount) || ...
             (isfinite(opts.magnitudeCoilCount) && ...
              opts.magnitudeCoilCount >= 1 && ...
              fix(opts.magnitudeCoilCount) == opts.magnitudeCoilCount)))
        error('T1dispersion:InvalidMagnitudeCoilCount', ...
            'magnitudeCoilCount must be NaN or a positive integer.');
    end
    polarityOptions = {'forcedNegativeIndices', 'forcedPositiveIndices'};
    for optionIndex = 1:numel(polarityOptions)
        optionName = polarityOptions{optionIndex};
        optionValue = opts.(optionName);
        if ~(isempty(optionValue) || isnumeric(optionValue) || iscell(optionValue))
            error('T1dispersion:InvalidPolarityConstraint', ...
                '%s must be numeric, a cell array, or empty.', optionName);
        end
    end
    if opts.rateGridPoints < 20 || opts.profileGridPoints < 50
        error('T1dispersion:GridTooSmall', ...
            'rateGridPoints must be >=20 and profileGridPoints >=50.');
    end
    if opts.confidenceLevel <= 0 || opts.confidenceLevel >= 1
        error('T1dispersion:InvalidConfidence', ...
            'confidenceLevel must lie strictly between 0 and 1.');
    end
    if ~isscalar(opts.boundCandidatePenaltyFraction) || ...
            ~isfinite(opts.boundCandidatePenaltyFraction) || ...
            opts.boundCandidatePenaltyFraction < 0
        error('T1dispersion:InvalidBoundPenalty', ...
            'boundCandidatePenaltyFraction must be a non-negative scalar.');
    end
    if ~isscalar(opts.maximumReliableSignTransitions) || ...
            ~isfinite(opts.maximumReliableSignTransitions) || ...
            opts.maximumReliableSignTransitions < 0 || ...
            fix(opts.maximumReliableSignTransitions) ~= opts.maximumReliableSignTransitions
        error('T1dispersion:InvalidSignTransitions', ...
            'maximumReliableSignTransitions must be a non-negative integer.');
    end
    if ~isscalar(opts.maximumSignIndexDisagreement) || ...
            ~isfinite(opts.maximumSignIndexDisagreement) || ...
            opts.maximumSignIndexDisagreement < 0 || ...
            fix(opts.maximumSignIndexDisagreement) ~= opts.maximumSignIndexDisagreement
        error('T1dispersion:InvalidSignIndexDifference', ...
            'maximumSignIndexDisagreement must be a non-negative integer.');
    end
    if ~(isscalar(opts.detectionField_mT) && ...
            (isnan(opts.detectionField_mT) || ...
             (isfinite(opts.detectionField_mT) && opts.detectionField_mT > 0)))
        error('T1dispersion:InvalidDetectionField', ...
            'detectionField_mT must be NaN or a positive finite scalar.');
    end
    integerJointOptions = {'jointMinimumFields', 'jointRateSweeps', ...
        'jointRobustIterations', 'jointSignIterations'};
    for optionIndex = 1:numel(integerJointOptions)
        optionName = integerJointOptions{optionIndex};
        optionValue = opts.(optionName);
        if ~isscalar(optionValue) || ~isfinite(optionValue) || ...
                optionValue < 1 || fix(optionValue) ~= optionValue
            error('T1dispersion:InvalidJointOption', ...
                '%s must be a positive integer.', optionName);
        end
    end
    if opts.jointMinimumFields < 2
        error('T1dispersion:InvalidJointFieldCount', ...
            'jointMinimumFields must be at least 2.');
    end
    if ~isscalar(opts.jointTolerance) || ~isfinite(opts.jointTolerance) || ...
            opts.jointTolerance <= 0
        error('T1dispersion:InvalidJointTolerance', ...
            'jointTolerance must be a positive finite scalar.');
    end
    if ~isscalar(opts.jointMinimumAlpha) || ...
            ~isfinite(opts.jointMinimumAlpha) || opts.jointMinimumAlpha < 0
        error('T1dispersion:InvalidJointAlphaFloor', ...
            'jointMinimumAlpha must be a non-negative finite scalar.');
    end
    if ~isscalar(opts.jointLateTimeWeight) || ...
            ~isfinite(opts.jointLateTimeWeight) || ...
            opts.jointLateTimeWeight < 1
        error('T1dispersion:InvalidLateTimeWeight', ...
            'jointLateTimeWeight must be a finite scalar >= 1.');
    end
    if ~isscalar(opts.jointTimeWeightExponent) || ...
            ~isfinite(opts.jointTimeWeightExponent) || ...
            opts.jointTimeWeightExponent <= 0
        error('T1dispersion:InvalidTimeWeightExponent', ...
            'jointTimeWeightExponent must be a positive finite scalar.');
    end
    if ~isscalar(opts.jointMaximumNoiseWeightRatio) || ...
            ~isfinite(opts.jointMaximumNoiseWeightRatio) || ...
            opts.jointMaximumNoiseWeightRatio < 1
        error('T1dispersion:InvalidJointNoiseRatio', ...
            'jointMaximumNoiseWeightRatio must be at least 1.');
    end
    logicalJointOptions = {'jointNoiseWeighting', ...
        'jointTimeWeighting', 'jointFallbackToIndependent'};
    for optionIndex = 1:numel(logicalJointOptions)
        optionName = logicalJointOptions{optionIndex};
        optionValue = opts.(optionName);
        if ~(islogical(optionValue) && isscalar(optionValue)) && ...
                ~(isnumeric(optionValue) && isscalar(optionValue) && ...
                  isfinite(optionValue) && any(optionValue == [0 1]))
            error('T1dispersion:InvalidJointOption', ...
                '%s must be a logical scalar.', optionName);
        end
        opts.(optionName) = logical(optionValue);
    end
end

function polarity = fieldPolarityConstraints(obj, opts, fieldIndex, sampleCount)
    negative = polarityOptionForField(opts.forcedNegativeIndices, fieldIndex);
    positive = polarityOptionForField(opts.forcedPositiveIndices, fieldIndex);
    if memberExists(obj, 'manual_polarity_flips') && ...
            ~isempty(obj.manual_polarity_flips)
        manualNegative = polarityOptionForField( ...
            obj.manual_polarity_flips, fieldIndex);
        negative = [negative(:); manualNegative(:)];
    end
    negative = unique(round(double(negative(:).')));
    positive = unique(round(double(positive(:).')));
    negative = negative(isfinite(negative) & negative >= 1 & ...
        negative <= sampleCount);
    positive = positive(isfinite(positive) & positive >= 1 & ...
        positive <= sampleCount);
    overlap = intersect(negative, positive);
    if ~isempty(overlap)
        error('T1dispersion:ConflictingPolarityConstraint', ...
            ['Field %d constrains sample(s) %s to be both negative ', ...
             'and positive.'], fieldIndex, mat2str(overlap));
    end
    polarity = struct('negativeIndices', negative, ...
        'positiveIndices', positive);
end

function indices = polarityOptionForField(value, fieldIndex)
    indices = [];
    if isempty(value)
        return
    end
    if iscell(value)
        if fieldIndex <= numel(value) && ~isempty(value{fieldIndex})
            indices = value{fieldIndex};
        end
    elseif isnumeric(value)
        indices = value;
    end
    if ~isnumeric(indices)
        error('T1dispersion:InvalidPolarityConstraint', ...
            'Polarity indices for field %d must be numeric.', fieldIndex);
    end
end

function mask2d = selectMask(mask, slice)
    if ndims(mask) >= 3 && size(mask, 3) > 1
        if slice > size(mask, 3)
            error('T1dispersion:MaskSliceOutOfRange', ...
                'Mask has %d slices; slice %d was requested.', size(mask, 3), slice);
        end
        mask2d = mask(:, :, slice);
    else
        mask2d = mask(:, :, 1);
    end
    mask2d = logical(mask2d);
end

function sources = resolveSignalSources(obj, requestedMode, nFields, mask, slice)
    sources = struct('mode', '', 'useComplex', false, 'useMagnitude', false, ...
        'complexData', [], 'magnitudeData', [], ...
        'complexMember', '', 'magnitudeMember', '', ...
        'complexRejectedReason', '');

    hasComplex = memberExists(obj, 'compleximage') && ~isempty(obj.compleximage);
    hasMagnitude = memberExists(obj, 'magimage') && ~isempty(obj.magimage);
    complexFull = false;
    magnitudeFull = false;

    % buildimages intentionally leaves compleximage six-dimensional in RSSQ
    % mode, while magimage is the fully combined five-dimensional dataset.
    % In auto mode an uncombined complex array must therefore be ignored, not
    % allowed to prevent the valid magnitude dispersion from being fitted.
    if hasMagnitude
        validateImageGeometry(obj.magimage, mask, slice, 'magimage');
        magnitudeFull = size(obj.magimage, 5) >= nFields;
        if magnitudeFull
            sources.magnitudeMember = 'magimage';
        end
    end
    if hasComplex
        if ndims(obj.compleximage) > 5 && size(obj.compleximage, 6) > 1
            sources.complexRejectedReason = sprintf( ...
                'compleximage retains %d receiver channels', ...
                size(obj.compleximage, 6));
            if strcmp(requestedMode, 'complex')
                error('T1dispersion:UncombinedComplexImage', ...
                    ['signalMode="complex" was requested, but obj.compleximage ', ...
                     'still has %d receiver channels. Use an Adaptive or ', ...
                     'ACS/Walsh complex combination, or use auto/magnitude.'], ...
                    size(obj.compleximage, 6));
            end
            warning('T1dispersion:UncombinedComplexImageIgnored', ...
                ['Ignoring %s and fitting the complete magimage dispersion.'], ...
                sources.complexRejectedReason);
        else
            validateImageGeometry(obj.compleximage, mask, slice, 'compleximage');
            complexFull = size(obj.compleximage, 5) >= nFields;
            if complexFull
                sources.complexMember = 'compleximage';
            elseif size(obj.compleximage, 5) == 1 && nFields > 1
                warning('T1dispersion:DisplayedComplexImageIgnored', ...
                    ['obj.compleximage contains only one field and therefore ', ...
                     'appears to be the displayed-field image. It will not be ', ...
                     'used as though it represented the full dispersion data.']);
            end
        end
    end

    switch requestedMode
        case 'complex'
            if ~complexFull
                error('T1dispersion:IncompleteComplexFields', ...
                    ['Complex fitting requires [X Y slice time field] data ', ...
                     'for all %d fieldpoints; obj.compleximage contains %d field(s).'], ...
                    nFields, imageFieldCount(obj, 'compleximage'));
            end
            sources.mode = 'complex';
            sources.useComplex = true;
            sources.complexData = obj.compleximage;

        case 'magnitude'
            if complexFull
                sources.mode = 'magnitude';
                sources.useComplex = true; % magnitude after complex ROI mean
                sources.complexData = obj.compleximage;
            elseif magnitudeFull
                sources.mode = 'magnitude';
                sources.useMagnitude = true;
                sources.magnitudeData = obj.magimage;
            else
                error('T1dispersion:IncompleteMagnitudeFields', ...
                    'No image member contains all %d magnetic fields.', nFields);
            end

        otherwise % auto: fit both when possible
            if complexFull
                sources.useComplex = true;
                sources.complexData = obj.compleximage;
            end
            if magnitudeFull && ~complexFull
                sources.useMagnitude = true;
                sources.magnitudeData = obj.magimage;
            end

            if sources.useComplex
                sources.mode = 'hybrid_auto';
            elseif sources.useMagnitude
                sources.mode = 'magnitude';
            else
                error('T1dispersion:NoFullDispersionImage', ...
                    ['Neither obj.compleximage nor obj.magimage contains all ', ...
                     '%d fields in obj.fieldpoints. The fitter will not ', ...
                     'silently process only the displayed field.'], nFields);
            end
    end
end

function validateImageGeometry(imageData, mask, slice, memberName)
    if ndims(imageData) < 4
        error('T1dispersion:ImageDimensions', ...
            'obj.%s must have at least [X Y slice time] dimensions.', memberName);
    end
    if ndims(imageData) > 5 && size(imageData, 6) > 1
        error('T1dispersion:UncombinedCoils', ...
            'obj.%s contains an uncombined sixth coil dimension.', memberName);
    end
    if slice < 1 || slice > size(imageData, 3)
        error('T1dispersion:SliceOutOfRange', ...
            'Requested slice %d; obj.%s contains %d slice(s).', ...
            slice, memberName, size(imageData, 3));
    end
    if size(mask, 1) ~= size(imageData, 1) || size(mask, 2) ~= size(imageData, 2)
        error('T1dispersion:MaskSizeMismatch', ...
            'Mask size %s does not match obj.%s size [%d %d].', ...
            mat2str(size(mask)), memberName, size(imageData, 1), size(imageData, 2));
    end
end

function count = imageFieldCount(obj, memberName)
    count = 0;
    if memberExists(obj, memberName) && ~isempty(obj.(memberName))
        count = size(obj.(memberName), 5);
    end
end

function dimensions = memberSize(obj, memberName)
    if memberExists(obj, memberName) && ~isempty(obj.(memberName))
        dimensions = size(obj.(memberName));
    else
        dimensions = [];
    end
end

function coilCount = inferMagnitudeCoilCount(obj, sources)
    coilCount = 1;
    if isfield(sources, 'complexData') && ~isempty(sources.complexData) && ...
            ndims(sources.complexData) > 5 && size(sources.complexData, 6) > 1
        coilCount = size(sources.complexData, 6);
        return
    end
    if memberExists(obj, 'compleximage') && ~isempty(obj.compleximage) && ...
            ndims(obj.compleximage) > 5 && size(obj.compleximage, 6) > 1
        coilCount = size(obj.compleximage, 6);
        return
    end
    possibleMembers = {'n_coils', 'ncoils', 'channels'};
    for memberIndex = 1:numel(possibleMembers)
        name = possibleMembers{memberIndex};
        if memberExists(obj, name)
            value = double(obj.(name));
            if isscalar(value) && isfinite(value) && value >= 1
                coilCount = round(value);
                return
            end
        end
    end
end

function t = timeVectorForField(times, fieldIndex, nFields)
    if isvector(times)
        t = times(:);
    elseif size(times, 1) == nFields
        t = times(fieldIndex, :).';
    elseif size(times, 2) == nFields
        t = times(:, fieldIndex);
    else
        error('T1dispersion:TimeDimensions', ...
            ['timepoints size %s has neither rows nor columns equal to the ', ...
             '%d values in fieldpoints.'], mat2str(size(times)), nFields);
    end
    t = double(t(:));
end

function [t_s, signal] = prepareFieldSamples(t_ms, signal, fieldIndex, label)
    if isempty(signal)
        t_s = [];
        return
    end
    nUse = min(numel(t_ms), numel(signal));
    if numel(t_ms) ~= numel(signal)
        warning('T1dispersion:TimeImageMismatch', ...
            ['Field %d %s data: %d evolution times but %d images. ', ...
             'Using the first %d paired samples.'], ...
            fieldIndex, label, numel(t_ms), numel(signal), nUse);
    end
    t_s = double(t_ms(1:nUse)) ./ 1000;
    signal = signal(1:nUse);
    [t_s, signal] = cleanAndCombineSamples(t_s, signal);
end

function [signal, meta] = extractRoiSignal(imageData, mask, slice, fieldIndex, mode, opts)
    nx = size(imageData, 1);
    ny = size(imageData, 2);
    nt = size(imageData, 4);
    block = reshape(imageData(:, :, slice, :, fieldIndex), nx, ny, nt);
    allVoxels = reshape(block, nx * ny, nt);
    backgroundVoxels = allVoxels(~mask(:), :);
    [backgroundSigma, backgroundMean, backgroundSource] = ...
        estimateBackgroundNoise(backgroundVoxels, mode, opts);
    voxels = allVoxels(mask(:), :);

    meta = struct('voxelCount', size(voxels, 1), ...
        'usedVoxelCount', size(voxels, 1), 'anchorTimeIndex', NaN, ...
        'aggregation', opts.roiAggregation, ...
        'backgroundNoiseSigma', backgroundSigma, ...
        'backgroundMagnitudeMean', backgroundMean, ...
        'roiNoiseSigma', backgroundSigma, ...
        'noiseSource', backgroundSource);

    if strcmp(mode, 'magnitude')
        signal = nan(nt, 1);
        for k = 1:nt
            values = real(voxels(:, k));
            values = values(isfinite(values));
            signal(k) = robustSpatialMean(values, opts.trimFraction, ...
                strcmp(opts.roiAggregation, 'phase_aligned_trimmed_mean'));
        end
        return
    end

    if strcmp(opts.roiAggregation, 'plain_mean')
        signal = nan(nt, 1);
        for k = 1:nt
            values = voxels(:, k);
            good = isfinite(real(values)) & isfinite(imag(values));
            values = values(good);
            if ~isempty(values)
                signal(k) = mean(values);
            end
        end
        return
    end

    anchorScore = nan(1, nt);
    for k = 1:nt
        a = abs(voxels(:, k));
        a = a(isfinite(a));
        if ~isempty(a)
            anchorScore(k) = median(a);
        end
    end
    if ~any(isfinite(anchorScore))
        signal = complex(nan(nt, 1), nan(nt, 1));
        return
    end
    [~, anchorIndex] = max(anchorScore);
    anchor = voxels(:, anchorIndex);
    anchorMagnitude = abs(anchor);
    finiteAnchor = isfinite(real(anchor)) & isfinite(imag(anchor)) & ...
        isfinite(anchorMagnitude) & anchorMagnitude > 0;
    positiveAnchor = anchorMagnitude(finiteAnchor);
    if isempty(positiveAnchor)
        signal = complex(nan(nt, 1), nan(nt, 1));
        return
    end
    anchorThreshold = opts.anchorMinFraction * median(positiveAnchor);
    keepVoxel = finiteAnchor & anchorMagnitude >= anchorThreshold;
    voxels = voxels(keepVoxel, :);
    anchor = anchor(keepVoxel);
    rotation = conj(anchor ./ abs(anchor));
    voxels = bsxfun(@times, voxels, rotation);

    signal = complex(nan(nt, 1), nan(nt, 1));
    useTrim = strcmp(opts.roiAggregation, 'phase_aligned_trimmed_mean');
    for k = 1:nt
        values = voxels(:, k);
        good = isfinite(real(values)) & isfinite(imag(values));
        values = values(good);
        signal(k) = robustComplexSpatialMean(values, opts.trimFraction, useTrim);
    end
    meta.usedVoxelCount = size(voxels, 1);
    meta.anchorTimeIndex = anchorIndex;
    if isfinite(backgroundSigma) && backgroundSigma > 0 && ...
            meta.usedVoxelCount > 0
        retainedFraction = max(1 - opts.trimFraction, eps);
        meta.roiNoiseSigma = backgroundSigma ./ ...
            sqrt(meta.usedVoxelCount .* retainedFraction);
    end
end

function [sigma, magnitudeMean, source] = estimateBackgroundNoise( ...
        backgroundVoxels, mode, opts)
    sigma = NaN;
    magnitudeMean = NaN;
    source = 'unavailable';
    if isempty(backgroundVoxels)
        return
    end
    temporalScore = nan(size(backgroundVoxels, 1), 1);
    for backgroundIndex = 1:size(backgroundVoxels, 1)
        values = backgroundVoxels(backgroundIndex, :);
        good = isfinite(real(values)) & isfinite(imag(values));
        if any(good)
            temporalScore(backgroundIndex) = median(abs(values(good)));
        end
    end
    validScore = temporalScore(isfinite(temporalScore));
    if isempty(validScore)
        return
    end
    sortedScore = sort(validScore, 'ascend');
    selectionIndex = max(1, min(numel(sortedScore), ...
        ceil(opts.backgroundFraction .* numel(sortedScore))));
    threshold = sortedScore(selectionIndex);
    selected = backgroundVoxels(temporalScore <= threshold, :);
    selected = selected(:);
    selected = selected(isfinite(real(selected)) & isfinite(imag(selected)));
    if isempty(selected)
        return
    end
    magnitudeMean = mean(abs(selected));
    if strcmp(mode, 'complex')
        realScale = robustRealScale(real(selected));
        imaginaryScale = robustRealScale(imag(selected));
        validScale = [realScale imaginaryScale];
        validScale = validScale(isfinite(validScale) & validScale > 0);
        if ~isempty(validScale)
            sigma = median(validScale);
            source = 'low-signal complex background';
        end
    else
        magnitudeMedian = median(abs(selected));
        if isfinite(magnitudeMedian) && magnitudeMedian > 0
            % Exact for single-coil Rayleigh background.  For a combined
            % magnitude image this is used only as a scale estimate; the
            % noncentral-chi approximation uses the measured mean floor.
            sigma = magnitudeMedian ./ sqrt(2 .* log(2));
            source = 'low-signal magnitude background';
        end
    end
end

function value = robustSpatialMean(values, trimFraction, useTrim)
    if isempty(values)
        value = NaN;
        return
    end
    values = values(:);
    if ~useTrim || trimFraction <= 0 || numel(values) < 8
        value = mean(values);
        return
    end
    centre = median(values);
    [~, order] = sort(abs(values - centre), 'ascend');
    nKeep = max(1, floor((1 - trimFraction) * numel(values)));
    value = mean(values(order(1:nKeep)));
end

function value = robustComplexSpatialMean(values, trimFraction, useTrim)
    if isempty(values)
        value = complex(NaN, NaN);
        return
    end
    values = values(:);
    if ~useTrim || trimFraction <= 0 || numel(values) < 8
        value = mean(values);
        return
    end
    centre = median(real(values)) + 1i * median(imag(values));
    [~, order] = sort(abs(values - centre), 'ascend');
    nKeep = max(1, floor((1 - trimFraction) * numel(values)));
    value = mean(values(order(1:nKeep)));
end

function [t, y] = cleanAndCombineSamples(t, y)
    t = double(t(:));
    y = y(:);
    good = isfinite(t) & t >= 0 & isfinite(real(y)) & isfinite(imag(y));
    t = t(good);
    y = y(good);
    [t, order] = sort(t, 'ascend');
    y = y(order);
    if isempty(t)
        return
    end

    [uniqueTimes, ~, group] = unique(t);
    if numel(uniqueTimes) == numel(t)
        return
    end
    combined = zeros(size(uniqueTimes), 'like', y);
    for k = 1:numel(uniqueTimes)
        combined(k) = mean(y(group == k));
    end
    t = uniqueTimes;
    y = combined;
end


% ========================================================================
% Complex and magnitude fitting
% ========================================================================

function r = fitOrFailComplex(t, signal, opts, template, polarity, fieldIndex)
    if numel(t) < opts.minimumPoints
        r = template;
        r.method = 'complex_phase_drift_ir';
        r.tSeconds = t;
        r.rawSignal = signal;
        r.displaySignal = real(signal);
        r.qc = sprintf('failed: fewer than %d valid complex time points', ...
            opts.minimumPoints);
        r.qcReasons = {r.qc};
    else
        r = fitComplexPhaseDriftIR(t, signal, opts, polarity, fieldIndex);
    end
end

function r = fitOrFailMagnitude(t, signal, opts, template, polarity, ...
        fieldIndex, magnitudeMeta)
    if numel(t) < opts.minimumPoints
        r = template;
        r.method = 'magnitude_rician_ir';
        r.tSeconds = t;
        r.rawSignal = signal;
        r.displaySignal = real(signal);
        r.qc = sprintf('failed: fewer than %d valid magnitude time points', ...
            opts.minimumPoints);
        r.qcReasons = {r.qc};
    else
        r = fitMagnitudeRicianIR(t, signal, opts, polarity, ...
            fieldIndex, magnitudeMeta);
    end
end

function selected = selectHybridFit(complexFit, magnitudeFit, opts)
    complexProblems = {};
    if ~complexFit.ok
        complexProblems{end + 1} = 'complex fit failed'; %#ok<AGROW>
    else
        if complexFit.boundHit
            complexProblems{end + 1} = 'complex rate reached a bound'; %#ok<AGROW>
        end
        if isfinite(complexFit.phaseResidualDegrees) && ...
                complexFit.phaseResidualDegrees > ...
                opts.maximumPhaseResidualDegrees
            complexProblems{end + 1} = sprintf( ...
                'complex phase residual is %.3g degrees', ...
                complexFit.phaseResidualDegrees); %#ok<AGROW>
        end
        if complexFit.reliableSignTransitions > opts.maximumReliableSignTransitions
            complexProblems{end + 1} = sprintf( ...
                'complex projection has %d reliable sign transitions', ...
                complexFit.reliableSignTransitions); %#ok<AGROW>
        end
        if isfinite(complexFit.orthogonalNRMSE) && ...
                complexFit.orthogonalNRMSE > opts.maximumOrthogonalNRMSE
            complexProblems{end + 1} = 'complex phase trajectory is incoherent'; %#ok<AGROW>
        end
        if isfinite(complexFit.nrmse) && complexFit.nrmse > opts.maximumNRMSE
            complexProblems{end + 1} = 'complex residual is large'; %#ok<AGROW>
        end
        magnitudeBranchTrusted = magnitudeFit.ok && ~magnitudeFit.boundHit && ...
            ~magnitudeFit.signAmbiguous && isfinite(magnitudeFit.nrmse) && ...
            magnitudeFit.nrmse <= opts.maximumNRMSE;
        if magnitudeBranchTrusted && isfinite(complexFit.signChangeIndex) && ...
                isfinite(magnitudeFit.signChangeIndex) && ...
                abs(complexFit.signChangeIndex - magnitudeFit.signChangeIndex) > ...
                opts.maximumSignIndexDisagreement
            complexProblems{end + 1} = sprintf( ...
                'complex and magnitude fits place the zero crossing at indices %d and %d', ...
                complexFit.signChangeIndex, magnitudeFit.signChangeIndex); %#ok<AGROW>
        end
    end

    if isempty(complexProblems)
        selected = complexFit;
        selected.selectionReason = ...
            'complex fit retained: phase projection has at most one reliable sign transition';
    elseif magnitudeFit.ok
        selected = magnitudeFit;
        selected.selectionReason = ['magnitude-constrained fit selected: ' ...
            strjoin(complexProblems, '; ')];
    else
        selected = complexFit;
        selected.selectionReason = ['complex fit retained because magnitude fit failed: ' ...
            strjoin(complexProblems, '; ')];
    end

    selected.complexCandidateSummary = fitSummary(complexFit);
    selected.magnitudeCandidateSummary = fitSummary(magnitudeFit);
end

function summary = fitSummary(r)
    summary = struct('ok', r.ok, 'method', r.method, 'R1', r.R1, ...
        'T1ms', r.T1ms, 'boundHit', r.boundHit, ...
        'signChangeIndex', r.signChangeIndex, ...
        'reliableSignTransitions', r.reliableSignTransitions, ...
        'signAmbiguous', r.signAmbiguous, 'nrmse', r.nrmse, ...
        'orthogonalNRMSE', r.orthogonalNRMSE, ...
        'phaseResidualDegrees', r.phaseResidualDegrees, 'qc', r.qc);
end

function r = fitComplexPhaseDriftIR(t, z, opts, polarity, fieldIndex)
    if strcmp(opts.sequenceMode, 'unconstrained')
        r = fitComplexAffine(t, z, opts);
        return
    end
    r = emptyFieldResult();
    r.method = 'complex_phase_drift_ir';
    r.tSeconds = t(:);
    r.rawSignal = z(:);
    n = numel(t);
    signalScale = max(abs(z));
    if ~isfinite(signalScale) || signalScale <= eps
        r.qc = 'failed: complex signal has no finite dynamic range';
        r.qcReasons = {r.qc};
        return
    end

    candidateTemplate = struct('ok', false, 'score', Inf, 'R1', NaN, ...
        'fit', struct(), 'corrected', [], 'phaseTrend', [], ...
        'phaseResidualDegrees', Inf, 'branchIndex', NaN);
    candidates = repmat(candidateTemplate, 1, n + 1);
    magnitude = abs(z(:));
    phaseNoiseProxy = max(min(magnitude) ./ sqrt(pi ./ 2), eps);
    reliableThreshold = max(opts.phaseReliableMagnitudeFraction .* ...
        signalScale, opts.minimumPhaseSNR .* phaseNoiseProxy);
    reliable = magnitude >= reliableThreshold & isfinite(magnitude);
    if sum(reliable) < max(3, opts.phaseDriftOrder + 1)
        [~, reliableOrder] = sort(magnitude, 'descend');
        reliable(:) = false;
        reliable(reliableOrder(1:min(n, max(3, ...
            opts.phaseDriftOrder + 1)))) = true;
    end

    for branchIndex = 0:n
        if ~phaseBranchAllowed(branchIndex, n, polarity)
            continue
        end
        signs = ones(n, 1);
        signs(1:branchIndex) = -1;
        aligned = z(:) .* signs;
        [phaseTrend, phaseResidualDegrees] = smoothPhaseTrend( ...
            aligned, reliable, magnitude, opts.phaseDriftOrder);
        corrected = z(:) .* exp(-1i .* phaseTrend);
        branchPolarity = struct('negativeIndices', 1:branchIndex, ...
            'positiveIndices', (branchIndex + 1):n);
        signedFit = fitConstrainedSignedIR(t, real(corrected), ...
            ones(n, 1), opts, branchPolarity);
        if ~signedFit.ok
            continue
        end
        complexResidual = corrected - signedFit.prediction;
        score = sum(abs(complexResidual).^2) ./ max(signalScale.^2, eps);
        candidates(branchIndex + 1).ok = true;
        candidates(branchIndex + 1).score = score;
        candidates(branchIndex + 1).R1 = signedFit.R1;
        candidates(branchIndex + 1).fit = signedFit;
        candidates(branchIndex + 1).corrected = corrected;
        candidates(branchIndex + 1).phaseTrend = phaseTrend;
        candidates(branchIndex + 1).phaseResidualDegrees = ...
            phaseResidualDegrees;
        candidates(branchIndex + 1).branchIndex = branchIndex;
    end

    scores = [candidates.score];
    [sortedScores, order] = sort(scores, 'ascend');
    if isempty(order) || ~isfinite(sortedScores(1))
        r.qc = sprintf(['failed: no phase-drift IR branch satisfied the ', ...
            'polarity constraints for field %d'], fieldIndex);
        r.qcReasons = {r.qc};
        return
    end
    best = candidates(order(1));
    fit = best.fit;
    corrected = best.corrected;
    prediction = fit.prediction;
    residual = corrected - prediction;
    noiseSigma = robustComplexScale(residual);
    dfe = max(2 .* n - (4 + opts.phaseDriftOrder), 1);
    sse = sum(abs(residual).^2);
    rmse = sqrt(sse ./ dfe);
    sst = sum(abs(corrected - mean(corrected)).^2);
    scaleAmplitude = max([fit.C, signalScale, eps]);
    signAmbiguous = false;
    if numel(order) > 1 && isfinite(sortedScores(2))
        second = candidates(order(2));
        closeScore = sortedScores(2) <= opts.signAmbiguityRatio .* ...
            max(sortedScores(1), eps);
        differentRate = isfinite(second.R1) && ...
            abs(log(second.R1 ./ best.R1)) > ...
            opts.signAmbiguityLogRateDifference;
        signAmbiguous = closeScore && differentRate;
    end

    r.ok = true;
    r.B = fit.B;
    r.C = fit.C;
    r.S0 = fit.B - fit.C;
    r.Sinf = fit.B;
    r.R1 = fit.R1;
    r.T1ms = 1000 ./ fit.R1;
    r.R1CI95 = fit.R1CI95;
    r.T1CI95 = rateCiToT1Ci(fit.R1CI95);
    zCritical = sqrt(2) .* erfinv(opts.confidenceLevel);
    r.R1SE = intervalToSE(r.R1CI95, zCritical);
    r.T1SE = intervalToSE(r.T1CI95, zCritical);
    r.ciOpen = fit.ciOpen;
    if any(r.ciOpen)
        r.R1SE = NaN;
        r.T1SE = NaN;
    end
    r.phaseRad = median(best.phaseTrend);
    r.phaseTrendRad = best.phaseTrend;
    r.phaseResidualDegrees = best.phaseResidualDegrees;
    r.displaySignal = real(corrected);
    r.displayFit = prediction;
    r.predictionComplex = prediction;
    r.residual = residual;
    r.weights = fit.weights;
    r.outlierIndices = find(fit.weights < opts.outlierWeightThreshold).';
    r.noiseSigma = noiseSigma;
    r.dynamicSNR = fit.C ./ max(noiseSigma, eps);
    r.sse = sse;
    r.rsquare = safeRsquare(sse, sst);
    r.rmse = rmse;
    r.nrmse = rmse ./ scaleAmplitude;
    r.dfe = dfe;
    r.boundHit = fit.boundHit;
    r.zeroCrossingSeconds = zeroCrossingTime(fit.B, fit.C, fit.R1);
    r.signChangeIndex = best.branchIndex;
    r.reliableSignTransitions = countReliableSignTransitions( ...
        real(corrected), noiseSigma, fit.C);
    r.orthogonalOffset = 0;
    r.orthogonalNRMSE = sqrt(mean(imag(corrected).^2)) ./ scaleAmplitude;
    r.signAmbiguous = signAmbiguous;
    r.candidateRates = [candidates.R1];
    r.candidateScores = scores;
    r.selectionReason = sprintf( ...
        ['phase drift order %d; branch %d; RMS phase residual %.3g degrees'], ...
        opts.phaseDriftOrder, best.branchIndex, best.phaseResidualDegrees);
    r = applyQualityControl(r, opts);
end

function allowed = phaseBranchAllowed(branchIndex, sampleCount, polarity)
    negative = polarity.negativeIndices;
    positive = polarity.positiveIndices;
    allowed = all(negative <= branchIndex) && ...
        all(positive > branchIndex) && branchIndex >= 0 && ...
        branchIndex <= sampleCount;
end

function [trend, residualDegrees] = smoothPhaseTrend( ...
        alignedSignal, reliable, magnitude, requestedOrder)
    sampleCount = numel(alignedSignal);
    coordinate = linspace(-1, 1, sampleCount).';
    reliableIndex = find(reliable(:));
    order = min(requestedOrder, max(numel(reliableIndex) - 1, 0));
    designAll = polynomialDesign(coordinate, order);
    phase = unwrap(angle(alignedSignal(reliableIndex)));
    design = designAll(reliableIndex, :);
    weights = magnitude(reliableIndex).^2;
    weights = weights ./ max(mean(weights), eps);
    weightedDesign = bsxfun(@times, design, sqrt(weights));
    weightedPhase = phase .* sqrt(weights);
    coefficient = weightedDesign \ weightedPhase;
    trend = designAll * coefficient;
    circularResidual = angle(alignedSignal(reliableIndex) .* ...
        exp(-1i .* trend(reliableIndex)));
    residualDegrees = 180 ./ pi .* sqrt(sum(weights .* ...
        circularResidual.^2) ./ max(sum(weights), eps));
end

function design = polynomialDesign(coordinate, order)
    design = ones(numel(coordinate), order + 1);
    for polynomialOrder = 1:order
        design(:, polynomialOrder + 1) = coordinate(:) .^ polynomialOrder;
    end
end

function fit = fitConstrainedSignedIR(t, y, initialWeights, opts, polarity)
    weights = initialWeights(:);
    previousRate = NaN;
    for robustIndex = 1:max(1, opts.robustIterations)
        fit = optimiseConstrainedIRRate(t, y, weights, opts, polarity);
        if ~fit.ok || ~opts.robust
            break
        end
        residual = y(:) - fit.prediction;
        sigma = robustRealScale(residual);
        proposedWeights = huberWeights(abs(residual), sigma, opts);
        newWeights = 0.5 .* weights + 0.5 .* proposedWeights;
        if isfinite(previousRate)
            rateChange = abs(log(fit.R1 ./ previousRate));
        else
            rateChange = Inf;
        end
        weightChange = max(abs(newWeights - weights));
        weights = newWeights;
        previousRate = fit.R1;
        if rateChange < opts.robustTolerance && ...
                weightChange < opts.robustTolerance
            break
        end
    end
    fit = optimiseConstrainedIRRate(t, y, weights, opts, polarity);
    fit.weights = weights;
    if fit.ok
        [fit.R1CI95, fit.ciOpen] = constrainedRateCI( ...
            t, y, weights, fit.R1, opts, polarity);
    else
        fit.R1CI95 = [NaN NaN];
        fit.ciOpen = [false false];
    end
end

function fit = optimiseConstrainedIRRate(t, y, weights, opts, polarity)
    logLower = log(opts.R1Bounds(1));
    logUpper = log(opts.R1Bounds(2));
    logGrid = linspace(logLower, logUpper, opts.rateGridPoints);
    score = inf(size(logGrid));
    for rateIndex = 1:numel(logGrid)
        [~, ~, ~, score(rateIndex)] = solveConstrainedIRAtRate( ...
            exp(logGrid(rateIndex)), t, y, weights, opts, polarity);
    end
    [~, bestIndex] = min(score);
    candidateLogRate = [logGrid(bestIndex), logLower, logUpper];
    if bestIndex > 1 && bestIndex < numel(logGrid)
        fminOptions = optimset('Display', 'off', 'TolX', 1e-8, ...
            'MaxIter', 100, 'MaxFunEvals', 220);
        objective = @(logRate) constrainedIRProfileObjective( ...
            logRate, t, y, weights, opts, polarity);
        refined = fminbnd(objective, logGrid(bestIndex - 1), ...
            logGrid(bestIndex + 1), fminOptions);
        candidateLogRate(end + 1) = refined; %#ok<AGROW>
    end
    candidateScore = inf(size(candidateLogRate));
    for candidateIndex = 1:numel(candidateLogRate)
        candidateScore(candidateIndex) = constrainedIRProfileObjective( ...
            candidateLogRate(candidateIndex), t, y, weights, opts, polarity);
    end
    [~, selected] = min(candidateScore);
    selectedLogRate = candidateLogRate(selected);
    rate = exp(selectedLogRate);
    [B, C, prediction, sse] = solveConstrainedIRAtRate( ...
        rate, t, y, weights, opts, polarity);
    constraintsSatisfied = all(prediction(polarity.negativeIndices) < 0) && ...
        all(prediction(polarity.positiveIndices) > 0);
    logSpan = logUpper - logLower;
    fit = struct('ok', isfinite(sse) && isfinite(B) && isfinite(C) && ...
        constraintsSatisfied, ...
        'R1', rate, 'B', B, 'C', C, 'prediction', prediction, ...
        'sse', sse, 'boundHit', selectedLogRate <= ...
        logLower + 0.01 .* logSpan || selectedLogRate >= ...
        logUpper - 0.01 .* logSpan, 'weights', weights, ...
        'R1CI95', [NaN NaN], 'ciOpen', [false false]);
end

function value = constrainedIRProfileObjective( ...
        logRate, t, y, weights, opts, polarity)
    [~, ~, ~, value] = solveConstrainedIRAtRate( ...
        exp(logRate), t, y, weights, opts, polarity);
    if ~isfinite(value)
        value = realmax('double');
    end
end

function [B, C, prediction, sse] = solveConstrainedIRAtRate( ...
        rate, t, y, weights, opts, polarity)
    exponential = exp(-rate .* t(:));
    signalScale = max(max(abs(y)), eps);
    minimumInitial = opts.minimumInitialMagnetizationFraction .* signalScale;
    design = [1 - exponential, -exponential];
    adjustedObservation = y(:) + minimumInitial .* exponential;
    squareRootWeight = sqrt(max(weights(:), 0));
    weightedDesign = bsxfun(@times, design, squareRootWeight);
    weightedObservation = adjustedObservation .* squareRootWeight;
    coefficient = nonnegativeCoordinateLeastSquares( ...
        weightedDesign, weightedObservation, opts.robustTolerance);
    if numel(coefficient) ~= 2 || any(~isfinite(coefficient))
        B = NaN;
        C = NaN;
        prediction = nan(size(t(:)));
        sse = Inf;
        return
    end
    B = coefficient(1);
    initialMagnitude = minimumInitial + coefficient(2);
    C = B + initialMagnitude;
    prediction = B - C .* exponential;
    residual = y(:) - prediction;
    sse = sum(weights(:) .* residual.^2);
    scalePenalty = max(sum(weights(:) .* y(:).^2), signalScale.^2);
    negativeViolation = polarity.negativeIndices( ...
        prediction(polarity.negativeIndices) >= 0);
    positiveViolation = polarity.positiveIndices( ...
        prediction(polarity.positiveIndices) <= 0);
    if ~isempty(negativeViolation) || ~isempty(positiveViolation)
        sse = sse + 1e6 .* scalePenalty .* ...
            (numel(negativeViolation) + numel(positiveViolation));
    end
end

function [rateCI, openBound] = constrainedRateCI( ...
        t, y, weights, rateHat, opts, polarity)
    logGrid = linspace(log(opts.R1Bounds(1)), ...
        log(opts.R1Bounds(2)), opts.profileGridPoints);
    logGrid = unique(sort([logGrid log(rateHat)]));
    profile = inf(size(logGrid));
    for profileIndex = 1:numel(logGrid)
        profile(profileIndex) = constrainedIRProfileObjective( ...
            logGrid(profileIndex), t, y, weights, opts, polarity);
    end
    [minimumSse, minimumIndex] = min(profile);
    dfe = max(numel(t) - (4 + opts.phaseDriftOrder), 1);
    fThreshold = fQuantileOneDof(opts.confidenceLevel, dfe);
    threshold = minimumSse .* (1 + fThreshold ./ dfe);
    [rateCI, openBound] = intervalFromProfile( ...
        logGrid, profile, minimumIndex, threshold);
end

function quantile = fQuantileOneDof(probability, denominatorDof)
    try
        betaValue = betaincinv(probability, 0.5, denominatorDof ./ 2);
        quantile = denominatorDof .* betaValue ./ max(1 - betaValue, eps);
    catch
        quantile = 2 .* erfinv(probability).^2;
    end
end

function [rateCI, openBound] = intervalFromProfile( ...
        logGrid, profile, minimumIndex, threshold)
    inside = profile <= threshold;
    leftIndex = minimumIndex;
    while leftIndex > 1 && inside(leftIndex - 1)
        leftIndex = leftIndex - 1;
    end
    rightIndex = minimumIndex;
    while rightIndex < numel(logGrid) && inside(rightIndex + 1)
        rightIndex = rightIndex + 1;
    end
    leftOpen = leftIndex == 1;
    rightOpen = rightIndex == numel(logGrid);
    if leftOpen
        leftLog = logGrid(1);
    else
        leftLog = interpolateThreshold(logGrid(leftIndex - 1), ...
            profile(leftIndex - 1), logGrid(leftIndex), ...
            profile(leftIndex), threshold);
    end
    if rightOpen
        rightLog = logGrid(end);
    else
        rightLog = interpolateThreshold(logGrid(rightIndex), ...
            profile(rightIndex), logGrid(rightIndex + 1), ...
            profile(rightIndex + 1), threshold);
    end
    rateCI = sort(exp([leftLog rightLog]));
    openBound = [leftOpen rightOpen];
end

function r = fitMagnitudeRicianIR(t, magnitudeSignal, opts, polarity, ...
        fieldIndex, magnitudeMeta)
    if strcmp(opts.sequenceMode, 'unconstrained')
        r = fitMagnitudeAffine(t, magnitudeSignal, opts);
        return
    end
    r = emptyFieldResult();
    r.method = 'magnitude_rician_ir';
    r.tSeconds = t(:);
    magnitudeSignal = abs(real(magnitudeSignal(:)));
    r.rawSignal = magnitudeSignal;
    [noiseSigma, noiseFloor, noiseModel, noiseSource, coilCount] = ...
        resolveMagnitudeNoise(magnitudeSignal, opts, fieldIndex, magnitudeMeta);
    weights = ones(size(t(:)));
    previousRate = NaN;
    for robustIndex = 1:max(1, opts.robustIterations)
        fit = optimiseMagnitudeIRRate(t, magnitudeSignal, weights, ...
            noiseSigma, noiseFloor, noiseModel, opts, polarity);
        if ~fit.ok || ~opts.robust
            break
        end
        residual = magnitudeSignal - fit.magnitudePrediction;
        residualScale = robustRealScale(residual);
        proposedWeights = huberWeights(abs(residual), residualScale, opts);
        newWeights = 0.5 .* weights + 0.5 .* proposedWeights;
        if isfinite(previousRate)
            rateChange = abs(log(fit.R1 ./ previousRate));
        else
            rateChange = Inf;
        end
        weightChange = max(abs(newWeights - weights));
        weights = newWeights;
        previousRate = fit.R1;
        if rateChange < opts.robustTolerance && ...
                weightChange < opts.robustTolerance
            break
        end
    end
    fit = optimiseMagnitudeIRRate(t, magnitudeSignal, weights, ...
        noiseSigma, noiseFloor, noiseModel, opts, polarity);
    if ~fit.ok
        r.qc = sprintf('failed: constrained magnitude IR fit failed at field %d', ...
            fieldIndex);
        r.qcReasons = {r.qc};
        return
    end

    latentSignal = fit.latentPrediction;
    if strcmp(noiseModel, 'rician')
        debiasedMagnitude = sqrt(max(magnitudeSignal.^2 - ...
            2 .* noiseSigma.^2, 0));
    else
        debiasedMagnitude = sqrt(max(magnitudeSignal.^2 - noiseFloor.^2, 0));
    end
    latentSigns = sign(latentSignal);
    latentSigns(latentSigns == 0) = 1;
    signedDisplay = latentSigns .* debiasedMagnitude;
    residual = magnitudeSignal - fit.magnitudePrediction;
    residualScale = robustRealScale(residual);
    dfe = max(numel(t) - 3, 1);
    sse = sum(residual.^2);
    rmse = sqrt(sse ./ dfe);
    scaleAmplitude = max([fit.C, max(magnitudeSignal), eps]);
    sst = sum((magnitudeSignal - mean(magnitudeSignal)).^2);
    zCritical = sqrt(2) .* erfinv(opts.confidenceLevel);

    r.ok = true;
    r.B = fit.B;
    r.C = fit.C;
    r.S0 = fit.B - fit.C;
    r.Sinf = fit.B;
    r.R1 = fit.R1;
    r.T1ms = 1000 ./ fit.R1;
    r.R1CI95 = fit.R1CI95;
    r.T1CI95 = rateCiToT1Ci(fit.R1CI95);
    r.R1SE = intervalToSE(r.R1CI95, zCritical);
    r.T1SE = intervalToSE(r.T1CI95, zCritical);
    r.ciOpen = fit.ciOpen;
    if any(r.ciOpen)
        r.R1SE = NaN;
        r.T1SE = NaN;
    end
    r.displaySignal = signedDisplay;
    r.displayFit = latentSignal;
    r.magnitudePrediction = fit.magnitudePrediction;
    r.magnitudeNoiseModel = noiseModel;
    r.magnitudeNoiseSigma = noiseSigma;
    r.magnitudeNoiseFloor = noiseFloor;
    r.magnitudeCoilCount = coilCount;
    r.noiseSource = noiseSource;
    r.residual = residual;
    r.weights = weights;
    r.outlierIndices = find(weights < opts.outlierWeightThreshold).';
    r.noiseSigma = noiseSigma;
    r.dynamicSNR = fit.C ./ max(noiseSigma, eps);
    r.sse = sse;
    r.rsquare = safeRsquare(sse, sst);
    r.rmse = rmse;
    r.nrmse = rmse ./ scaleAmplitude;
    r.dfe = dfe;
    r.boundHit = fit.boundHit;
    r.zeroCrossingSeconds = zeroCrossingTime(fit.B, fit.C, fit.R1);
    r.signChangeIndex = sum(latentSignal < 0);
    r.reliableSignTransitions = double(r.signChangeIndex > 0 && ...
        r.signChangeIndex < numel(t));
    r.signAmbiguous = fit.signAmbiguous;
    r.candidateRates = fit.profileRates;
    r.candidateScores = fit.profileScores;
    r.selectionReason = sprintf( ...
        '%s magnitude IR; noise %s; sigma %.4g; %d effective coil(s)', ...
        noiseModel, noiseSource, noiseSigma, coilCount);
    r = applyQualityControl(r, opts);
end

function [sigma, noiseFloor, model, source, coilCount] = resolveMagnitudeNoise( ...
        signal, opts, fieldIndex, meta)
    sigma = numericOptionForField(opts.magnitudeNoiseSigma, fieldIndex);
    source = 'OPTIONS.magnitudeNoiseSigma';
    if ~isfinite(sigma) || sigma <= 0
        sigma = NaN;
        if isstruct(meta)
            if isfield(meta, 'aggregation') && ...
                    strcmp(meta.aggregation, 'magnitude_of_complex_roi_mean') && ...
                    isfield(meta, 'roiNoiseSigma') && ...
                    isfinite(meta.roiNoiseSigma) && meta.roiNoiseSigma > 0
                sigma = meta.roiNoiseSigma;
                source = 'complex background propagated through ROI mean';
            elseif isfield(meta, 'backgroundNoiseSigma') && ...
                    isfinite(meta.backgroundNoiseSigma) && ...
                    meta.backgroundNoiseSigma > 0
                sigma = meta.backgroundNoiseSigma;
                source = meta.noiseSource;
            end
        end
    end
    if ~isfinite(sigma) || sigma <= 0
        sigma = max(min(signal) ./ sqrt(pi ./ 2), ...
            0.01 .* max(signal));
        source = 'recovery minimum fallback';
    end

    coilCount = opts.magnitudeCoilCount;
    if isnan(coilCount)
        coilCount = 1;
    end
    model = opts.magnitudeNoiseModel;
    if strcmp(model, 'auto')
        if coilCount <= 1 || (isstruct(meta) && ...
                isfield(meta, 'aggregation') && ...
                strcmp(meta.aggregation, 'magnitude_of_complex_roi_mean'))
            model = 'rician';
            coilCount = 1;
        else
            model = 'noncentral_chi_approx';
        end
    end
    noiseFloor = NaN;
    if isstruct(meta) && isfield(meta, 'backgroundMagnitudeMean') && ...
            isfinite(meta.backgroundMagnitudeMean) && ...
            meta.backgroundMagnitudeMean > 0
        noiseFloor = meta.backgroundMagnitudeMean;
    end
    if ~isfinite(noiseFloor) || noiseFloor <= 0
        if strcmp(model, 'rician')
            noiseFloor = sigma .* sqrt(pi ./ 2);
        else
            noiseFloor = sigma .* sqrt(2 .* coilCount);
        end
    end
end

function value = numericOptionForField(option, fieldIndex)
    value = NaN;
    if isempty(option)
        return
    end
    if isscalar(option)
        value = double(option);
    elseif fieldIndex <= numel(option)
        value = double(option(fieldIndex));
    end
end

function fit = optimiseMagnitudeIRRate(t, magnitude, weights, sigma, ...
        noiseFloor, model, opts, polarity)
    logLower = log(opts.R1Bounds(1));
    logUpper = log(opts.R1Bounds(2));
    logGrid = linspace(logLower, logUpper, opts.rateGridPoints);
    profile = inf(size(logGrid));
    for profileIndex = 1:numel(logGrid)
        solution = solveMagnitudeIRAtRate(exp(logGrid(profileIndex)), ...
            t, magnitude, weights, sigma, noiseFloor, model, opts, polarity);
        profile(profileIndex) = solution.score;
    end
    [~, bestIndex] = min(profile);
    candidateLogRate = [logGrid(bestIndex), logLower, logUpper];
    if bestIndex > 1 && bestIndex < numel(logGrid)
        fminOptions = optimset('Display', 'off', 'TolX', 2e-6, ...
            'MaxIter', 80, 'MaxFunEvals', 180);
        objective = @(logRate) magnitudeIRProfileObjective(logRate, ...
            t, magnitude, weights, sigma, noiseFloor, model, opts, polarity);
        refined = fminbnd(objective, logGrid(bestIndex - 1), ...
            logGrid(bestIndex + 1), fminOptions);
        candidateLogRate(end + 1) = refined; %#ok<AGROW>
    end
    candidateScore = inf(size(candidateLogRate));
    for candidateIndex = 1:numel(candidateLogRate)
        candidateScore(candidateIndex) = magnitudeIRProfileObjective( ...
            candidateLogRate(candidateIndex), t, magnitude, weights, ...
            sigma, noiseFloor, model, opts, polarity);
    end
    [~, selected] = min(candidateScore);
    selectedLogRate = candidateLogRate(selected);
    selectedRate = exp(selectedLogRate);
    solution = solveMagnitudeIRAtRate(selectedRate, t, magnitude, ...
        weights, sigma, noiseFloor, model, opts, polarity);
    profileRates = exp(logGrid);
    [profileRates, profileOrder] = sort([profileRates selectedRate]);
    profileScores = [profile solution.score];
    profileScores = profileScores(profileOrder);
    [profileRates, uniqueIndex] = unique(profileRates, 'stable');
    profileScores = profileScores(uniqueIndex);
    threshold = solution.score + erfinv(opts.confidenceLevel).^2;
    [~, minimumIndex] = min(profileScores);
    [rateCI, ciOpen] = intervalFromProfile(log(profileRates), ...
        profileScores, minimumIndex, threshold);

    localMinimum = false(size(profileScores));
    if numel(profileScores) >= 3
        localMinimum(2:end-1) = profileScores(2:end-1) <= ...
            profileScores(1:end-2) & profileScores(2:end-1) <= ...
            profileScores(3:end);
    end
    alternative = find(localMinimum & profileScores <= threshold & ...
        abs(log(profileRates ./ selectedRate)) > ...
        opts.signAmbiguityLogRateDifference, 1);
    signAmbiguous = ~isempty(alternative);
    logSpan = logUpper - logLower;
    fit = solution;
    fit.ok = isfinite(solution.score) && isfinite(solution.B) && ...
        isfinite(solution.C) && solution.C > solution.B && ...
        solution.constraintsSatisfied;
    fit.R1 = selectedRate;
    fit.boundHit = selectedLogRate <= logLower + 0.01 .* logSpan || ...
        selectedLogRate >= logUpper - 0.01 .* logSpan;
    fit.R1CI95 = rateCI;
    fit.ciOpen = ciOpen;
    fit.signAmbiguous = signAmbiguous;
    fit.profileRates = profileRates;
    fit.profileScores = profileScores;
end

function value = magnitudeIRProfileObjective(logRate, t, magnitude, ...
        weights, sigma, noiseFloor, model, opts, polarity)
    solution = solveMagnitudeIRAtRate(exp(logRate), t, magnitude, ...
        weights, sigma, noiseFloor, model, opts, polarity);
    value = solution.score;
    if ~isfinite(value)
        value = realmax('double');
    end
end

function solution = solveMagnitudeIRAtRate(rate, t, magnitude, weights, ...
        sigma, noiseFloor, model, opts, polarity)
    scale = max(max(magnitude), eps);
    minimumInitial = opts.minimumInitialMagnetizationFraction .* scale;
    lateCount = min(3, numel(magnitude));
    B0 = max(mean(magnitude(end-lateCount+1:end)), 0.1 .* scale);
    initialStarts = [B0, max(magnitude(1), B0); ...
        0.5 .* scale, scale; scale, scale];
    bestScore = Inf;
    bestParameters = [NaN NaN];
    fminOptions = optimset('Display', 'off', 'TolX', 1e-7, ...
        'TolFun', 1e-8, 'MaxIter', 180, 'MaxFunEvals', 360);
    for startIndex = 1:size(initialStarts, 1)
        startB = max(initialStarts(startIndex, 1), 1e-6 .* scale);
        startAExcess = max(initialStarts(startIndex, 2) - ...
            minimumInitial, 1e-6 .* scale);
        start = log([startB startAExcess] ./ scale);
        objective = @(parameters) magnitudeAmplitudeObjective(parameters, ...
            rate, t, magnitude, weights, sigma, noiseFloor, model, ...
            opts, polarity, scale, minimumInitial);
        [parameters, score] = fminsearch(objective, start, fminOptions);
        if score < bestScore
            bestScore = score;
            bestParameters = parameters;
        end
    end
    [B, initialMagnitude] = decodeMagnitudeAmplitudes( ...
        bestParameters, scale, minimumInitial);
    C = B + initialMagnitude;
    latent = B - C .* exp(-rate .* t(:));
    predictedMagnitude = magnitudeObservationMean(abs(latent), sigma, ...
        noiseFloor, model);
    constraintsSatisfied = all(latent(polarity.negativeIndices) < 0) && ...
        all(latent(polarity.positiveIndices) > 0);
    solution = struct('score', bestScore, 'B', B, 'C', C, ...
        'latentPrediction', latent, ...
        'magnitudePrediction', predictedMagnitude, ...
        'constraintsSatisfied', constraintsSatisfied);
end

function score = magnitudeAmplitudeObjective(parameters, rate, t, magnitude, ...
        weights, sigma, noiseFloor, model, opts, polarity, scale, minimumInitial)
    [B, initialMagnitude] = decodeMagnitudeAmplitudes( ...
        parameters, scale, minimumInitial);
    C = B + initialMagnitude;
    latent = B - C .* exp(-rate .* t(:));
    negativeValues = latent(polarity.negativeIndices);
    positiveValues = latent(polarity.positiveIndices);
    constraintPenalty = 1e6 .* (sum(max(negativeValues, 0).^2) + ...
        sum(max(-positiveValues, 0).^2)) ./ max(scale.^2, eps);
    noncentrality = abs(latent);
    if strcmp(model, 'rician')
        argument = magnitude(:) .* noncentrality ./ max(sigma.^2, eps);
        logBessel = log(max(besseli(0, argument, 1), realmin)) + ...
            abs(argument);
        term = (magnitude(:).^2 + noncentrality.^2) ./ ...
            max(2 .* sigma.^2, eps) - logBessel;
        score = sum(weights(:) .* term);
    else
        predicted = magnitudeObservationMean(noncentrality, sigma, ...
            noiseFloor, model);
        varianceScale = max(sigma.^2, (0.02 .* scale).^2);
        score = 0.5 .* sum(weights(:) .* ...
            (magnitude(:) - predicted).^2 ./ varianceScale);
    end
    score = score + constraintPenalty;
    if ~isfinite(score)
        score = realmax('double');
    end
end

function [B, initialMagnitude] = decodeMagnitudeAmplitudes( ...
        parameters, scale, minimumInitial)
    parameters = min(max(real(parameters(:).'), -20), log(50));
    B = scale .* exp(parameters(1));
    initialMagnitude = minimumInitial + scale .* exp(parameters(2));
end

function meanMagnitude = magnitudeObservationMean( ...
        noncentrality, sigma, noiseFloor, model)
    if strcmp(model, 'rician')
        if ~isfinite(sigma) || sigma <= eps
            meanMagnitude = noncentrality;
            return
        end
        q = noncentrality.^2 ./ (4 .* sigma.^2);
        scaledI0 = besseli(0, q, 1);
        scaledI1 = besseli(1, q, 1);
        meanMagnitude = sigma .* sqrt(pi ./ 2) .* ...
            ((1 + 2 .* q) .* scaledI0 + 2 .* q .* scaledI1);
        large = q > 1e4 | ~isfinite(meanMagnitude);
        meanMagnitude(large) = noncentrality(large) + ...
            sigma.^2 ./ max(2 .* noncentrality(large), eps);
    else
        meanMagnitude = sqrt(noncentrality.^2 + noiseFloor.^2);
    end
end

function r = fitComplexAffine(t, z, opts)
    r = emptyFieldResult();
    r.method = 'complex_affine_variable_projection';
    r.tSeconds = t(:);
    r.rawSignal = z(:);

    weights = ones(size(t));
    previousR = NaN;
    for iteration = 1:max(1, opts.robustIterations)
        [R1, beta, ~, boundHit] = optimiseRate(t, z, weights, opts);
        prediction = affinePrediction(t, R1, beta);
        residual = z - prediction;
        noiseSigma = robustComplexScale(residual);

        if ~opts.robust
            break
        end
        newWeights = huberWeights(abs(residual), noiseSigma, opts);
        if isfinite(previousR)
            rateChange = abs(log(R1 / previousR));
        else
            rateChange = Inf;
        end
        weightChange = max(abs(newWeights - weights));
        weights = 0.5 * weights + 0.5 * newWeights;
        previousR = R1;
        if rateChange < opts.robustTolerance && weightChange < opts.robustTolerance
            break
        end
    end
    [R1, beta, ~, boundHit] = optimiseRate(t, z, weights, opts);
    prediction = affinePrediction(t, R1, beta);
    residual = z - prediction;
    noiseSigma = robustComplexScale(residual);

    % Choose the arbitrary global phase so the exponential coefficient is
    % negative.  The displayed scalar curve is then B - C exp(-R1 t).
    if abs(beta(2)) > 0
        phase = angle(-beta(2));
    elseif abs(beta(1)) > 0
        phase = angle(beta(1));
    else
        phase = 0;
    end
    rotation = exp(-1i * phase);
    rotatedBeta = beta .* rotation;
    B = real(rotatedBeta(1));
    C = max(-real(rotatedBeta(2)), 0);
    displaySignal = real(z .* rotation);
    displayFit = B - C .* exp(-R1 .* t);
    orthogonalOffset = imag(rotatedBeta(1));
    orthogonalResidual = imag(z .* rotation) - orthogonalOffset;

    n = numel(t);
    dfe = max(2 * n - 5, 1);
    sse = sum(abs(residual).^2);
    sst = sum(abs(z - mean(z)).^2);
    rmse = sqrt(sse / dfe);
    rsq = safeRsquare(sse, sst);
    scaleAmplitude = max([C, max(abs(displaySignal - median(displaySignal))), eps]);

    [R1CI, ciOpen] = profileRateCI(t, z, weights, R1, 2 * n, 5, opts);
    T1CI = rateCiToT1Ci(R1CI);
    zCritical = sqrt(2) * erfinv(opts.confidenceLevel);

    r.ok = isfinite(R1) && R1 > 0 && all(isfinite(beta));
    r.B = B;
    r.C = C;
    r.S0 = B - C;
    r.Sinf = B;
    r.R1 = R1;
    r.T1ms = 1000 / R1;
    r.R1CI95 = R1CI;
    r.T1CI95 = T1CI;
    r.R1SE = intervalToSE(R1CI, zCritical);
    r.T1SE = intervalToSE(T1CI, zCritical);
    if any(ciOpen)
        r.R1SE = NaN;
        r.T1SE = NaN;
    end
    r.betaComplex = beta(:).';
    r.phaseRad = phase;
    r.signChangeIndex = sum(displayFit < 0);
    r.reliableSignTransitions = countReliableSignTransitions( ...
        displaySignal, noiseSigma, C);
    r.orthogonalOffset = orthogonalOffset;
    r.orthogonalNRMSE = sqrt(mean(orthogonalResidual.^2)) / scaleAmplitude;
    r.predictionComplex = prediction;
    r.displaySignal = displaySignal;
    r.displayFit = displayFit;
    r.residual = residual;
    r.weights = weights;
    r.outlierIndices = find(weights < opts.outlierWeightThreshold).';
    r.noiseSigma = noiseSigma;
    r.dynamicSNR = C / max(noiseSigma, eps);
    r.sse = sse;
    r.rsquare = rsq;
    r.rmse = rmse;
    r.nrmse = rmse / scaleAmplitude;
    r.dfe = dfe;
    r.boundHit = boundHit;
    r.ciOpen = ciOpen;
    r.zeroCrossingSeconds = zeroCrossingTime(B, C, R1);
    r = applyQualityControl(r, opts);
end

function r = fitMagnitudeAffine(t, magnitudeSignal, opts)
    r = emptyFieldResult();
    r.method = 'magnitude_all_zero_crossings';
    r.tSeconds = t(:);
    magnitudeSignal = abs(real(magnitudeSignal(:)));
    r.rawSignal = magnitudeSignal;

    n = numel(t);
    candidate = repmat(emptyMagnitudeCandidate(), 1, n);
    for k = 0:(n - 1)
        signs = ones(n, 1);
        if k > 0
            signs(1:k) = -1;
        end
        candidate(k + 1) = fitSignedCandidate(t, signs .* magnitudeSignal, ...
            magnitudeSignal, signs, k, opts);
    end

    scores = [candidate.score];
    [sortedScores, order] = sort(scores, 'ascend');
    if isempty(order) || ~isfinite(sortedScores(1))
        r.qc = 'failed: no magnitude sign pattern produced a finite fit';
        r.qcReasons = {r.qc};
        return
    end
    best = candidate(order(1));

    signAmbiguous = false;
    if numel(order) >= 2 && isfinite(sortedScores(2))
        signalScale = max(sum(magnitudeSignal.^2), eps);
        closeScore = sortedScores(2) <= opts.signAmbiguityRatio * ...
            max(sortedScores(1), eps * signalScale);
        differentRate = abs(log(candidate(order(2)).R1 / best.R1)) > ...
            opts.signAmbiguityLogRateDifference;
        signAmbiguous = closeScore && differentRate;
    end

    beta = best.beta;
    ySigned = best.ySigned;
    prediction = best.prediction;
    % Overall polarity is unknowable from magnitude. Orient the curve so the
    % exponential coefficient is negative, matching B - C exp(-R1 t).
    if beta(2) > 0
        beta = -beta;
        ySigned = -ySigned;
        prediction = -prediction;
    end
    B = beta(1);
    C = max(-beta(2), 0);
    residual = ySigned - prediction;
    n = numel(t);
    dfe = max(n - 3, 1);
    sse = sum(residual.^2);
    sst = sum((ySigned - mean(ySigned)).^2);
    rmse = sqrt(sse / dfe);
    scaleAmplitude = max([C, max(abs(ySigned - median(ySigned))), eps]);
    noiseSigma = robustRealScale(residual);

    [R1CI, ciOpen] = profileRateCI(t, ySigned, best.weights, ...
        best.R1, n, 3, opts);
    T1CI = rateCiToT1Ci(R1CI);
    zCritical = sqrt(2) * erfinv(opts.confidenceLevel);

    r.ok = isfinite(best.R1) && best.R1 > 0;
    r.B = B;
    r.C = C;
    r.S0 = B - C;
    r.Sinf = B;
    r.R1 = best.R1;
    r.T1ms = 1000 / best.R1;
    r.R1CI95 = R1CI;
    r.T1CI95 = T1CI;
    r.R1SE = intervalToSE(R1CI, zCritical);
    r.T1SE = intervalToSE(T1CI, zCritical);
    if any(ciOpen)
        r.R1SE = NaN;
        r.T1SE = NaN;
    end
    r.displaySignal = ySigned;
    r.displayFit = prediction;
    r.residual = residual;
    r.weights = best.weights;
    r.outlierIndices = find(best.weights < opts.outlierWeightThreshold).';
    r.noiseSigma = noiseSigma;
    r.dynamicSNR = C / max(noiseSigma, eps);
    r.sse = sse;
    r.rsquare = safeRsquare(sse, sst);
    r.rmse = rmse;
    r.nrmse = rmse / scaleAmplitude;
    r.dfe = dfe;
    r.boundHit = best.boundHit;
    r.ciOpen = ciOpen;
    % signChangeIndex is the number of negative samples in the final,
    % increasing orientation. It is n for a wholly negative -M0-to-zero
    % recovery, zero for a wholly positive saturation recovery, and an
    % interior index for a genuine zero crossing.
    r.signChangeIndex = sum(ySigned < 0);
    if r.signChangeIndex > 0 && r.signChangeIndex < n
        r.reliableSignTransitions = 1;
    else
        r.reliableSignTransitions = 0;
    end
    r.signAmbiguous = signAmbiguous;
    r.zeroCrossingSeconds = zeroCrossingTime(B, C, best.R1);
    r.candidateRates = [candidate.R1];
    r.candidateScores = scores;
    r = applyQualityControl(r, opts);
end

function c = fitSignedCandidate(t, y, magnitudeSignal, desiredSigns, signIndex, opts)
    c = emptyMagnitudeCandidate();
    c.signChangeIndex = signIndex;
    weights = ones(size(t));
    previousR = NaN;
    for iteration = 1:max(1, opts.robustIterations)
        [R1, beta, ~, boundHit] = optimiseRate(t, y, weights, opts);
        prediction = affinePrediction(t, R1, beta);
        residual = y - prediction;
        noiseSigma = robustRealScale(residual);
        if ~opts.robust
            break
        end
        newWeights = huberWeights(abs(residual), noiseSigma, opts);
        if isfinite(previousR)
            rateChange = abs(log(R1 / previousR));
        else
            rateChange = Inf;
        end
        weightChange = max(abs(newWeights - weights));
        weights = 0.5 * weights + 0.5 * newWeights;
        previousR = R1;
        if rateChange < opts.robustTolerance && weightChange < opts.robustTolerance
            break
        end
    end
    [R1, beta, ~, boundHit] = optimiseRate(t, y, weights, opts);
    prediction = affinePrediction(t, R1, beta);
    residual = y - prediction;
    noiseSigma = robustRealScale(residual);
    % Use the ordinary residual energy for choosing the polarity branch.
    % Allowing every branch to hide different samples with robust weights can
    % otherwise reward an implausible sign pattern. Robust weights are still
    % retained for the final parameter estimate after the branch is selected.
    score = sum(residual.^2);

    % A candidate sign pattern is only meaningful if the fitted signed curve
    % follows it away from the noise floor. Penalise clear inconsistencies.
    strong = magnitudeSignal > 2 * max(noiseSigma, eps);
    predictedSigns = sign(prediction);
    predictedSigns(predictedSigns == 0) = desiredSigns(predictedSigns == 0);
    mismatch = strong & predictedSigns ~= desiredSigns;
    if any(mismatch)
        score = score + 4 * sum(magnitudeSignal(mismatch).^2);
    end
    if boundHit
        score = score + opts.boundCandidatePenaltyFraction * ...
            max(sum(magnitudeSignal.^2), eps);
    end
    if signIndex > 0 && beta(2) >= 0
        % An interior negative-to-positive branch must be increasing. A
        % non-negative exponential coefficient describes a decreasing curve
        % and is therefore not a physically admissible orientation.
        score = Inf;
    end

    c.R1 = R1;
    c.beta = beta;
    c.ySigned = y;
    c.prediction = prediction;
    c.weights = weights;
    c.score = score;
    c.boundHit = boundHit;
end

function count = countReliableSignTransitions(signal, noiseSigma, amplitude)
    signal = real(signal(:));
    threshold = max(2 * max(noiseSigma, eps), 0.05 * max(amplitude, eps));
    reliable = abs(signal) > threshold & isfinite(signal);
    signs = sign(signal(reliable));
    if numel(signs) < 2
        count = 0;
    else
        count = sum(signs(2:end) ~= signs(1:end-1));
    end
end


% ========================================================================
% Joint multi-field fitting
% ========================================================================

function diagnostics = emptyJointDiagnostics(evolutionFields_mT, opts)
    diagnostics = struct( ...
        'attempted', false, ...
        'used', false, ...
        'modelName', ...
            'joint_multifield_inversion_constrained_variable_projection', ...
        'modelEquation', ['S_f(t) = M0_D*[r_f*(1-exp(-R1_f*t)) - ', ...
            'alpha_f*exp(-R1_f*t)]'], ...
        'representation', 'per-field phase-projected/signed ROI signal', ...
        'paperRelation', ['M0_D=C*B_D; r_f=B_E,f/B_D; ', ...
            'S_f(0)=-M0_D*alpha_f'], ...
        'fields_mT', evolutionFields_mT(:).', ...
        'fieldRatios', nan(size(evolutionFields_mT(:).')), ...
        'detectionField_mT', opts.detectionField_mT, ...
        'detectionFieldInferred', isnan(opts.detectionField_mT), ...
        'minimumAlpha', opts.jointMinimumAlpha, ...
        'timeWeightingEnabled', opts.jointTimeWeighting, ...
        'configuredLateTimeWeight', opts.jointLateTimeWeight, ...
        'lateTimeWeight', 1 + double(opts.jointTimeWeighting) * ...
            (opts.jointLateTimeWeight - 1), ...
        'timeWeightExponent', opts.jointTimeWeightExponent, ...
        'timeWeights', {cell(1, numel(evolutionFields_mT))}, ...
        'eligibleFieldIndices', [], ...
        'excludedFieldIndices', [], ...
        'initialR1', nan(size(evolutionFields_mT(:).')), ...
        'finalR1', nan(size(evolutionFields_mT(:).')), ...
        'initialSharedM0Detection', NaN, ...
        'sharedM0Detection', NaN, ...
        'sharedCPer_mT', NaN, ...
        'inversionAmplitude', nan(size(evolutionFields_mT(:).')), ...
        'Q', nan(size(evolutionFields_mT(:).')), ...
        'alpha', nan(size(evolutionFields_mT(:).')), ...
        'noisePrecision', nan(size(evolutionFields_mT(:).')), ...
        'magnitudeBranchIndex', nan(size(evolutionFields_mT(:).')), ...
        'weightedSSE', NaN, ...
        'globalDfe', NaN, ...
        'normalMatrixRcond', NaN, ...
        'rateUncertaintyMethod', ...
            'approximate joint weighted-Jacobian covariance', ...
        'robustIterations', 0, ...
        'rateSweeps', 0, ...
        'optimizerFunctionEvaluations', 0, ...
        'converged', false, ...
        'failureIdentifier', '', ...
        'failureMessage', '');
end

function [updatedResults, diagnostics] = fitJointMultiField( ...
        initialResults, evolutionFields_mT, opts)
    diagnostics = emptyJointDiagnostics(evolutionFields_mT, opts);
    diagnostics.attempted = true;
    updatedResults = initialResults;

    evolutionFields_mT = double(evolutionFields_mT(:).');
    if any(~isfinite(evolutionFields_mT)) || any(evolutionFields_mT < 0)
        error('T1dispersion:InvalidJointFields', ...
            'Joint fitting requires finite, non-negative evolution fields.');
    end
    positiveFields = evolutionFields_mT(evolutionFields_mT > 0);
    if isempty(positiveFields)
        error('T1dispersion:NoPositiveJointField', ...
            'At least one positive evolution field is required.');
    end
    if isnan(opts.detectionField_mT)
        detectionField_mT = max(positiveFields);
    else
        detectionField_mT = opts.detectionField_mT;
    end
    fieldRatios = evolutionFields_mT ./ detectionField_mT;
    diagnostics.detectionField_mT = detectionField_mT;
    diagnostics.fieldRatios = fieldRatios;

    eligible = false(1, numel(initialResults));
    for jointEligibilityIndex = 1:numel(initialResults)
        candidateResult = initialResults(jointEligibilityIndex);
        eligible(jointEligibilityIndex) = ...
            numel(candidateResult.tSeconds) >= opts.minimumPoints && ...
            numel(candidateResult.displaySignal) == ...
                numel(candidateResult.tSeconds) && ...
            all(isfinite(candidateResult.tSeconds(:))) && ...
            all(isfinite(real(candidateResult.displaySignal(:))));
    end
    eligibleIndices = find(eligible);
    diagnostics.eligibleFieldIndices = eligibleIndices;
    diagnostics.excludedFieldIndices = find(~eligible);
    if numel(eligibleIndices) < opts.jointMinimumFields
        diagnostics.failureIdentifier = 'T1dispersion:InsufficientJointFields';
        diagnostics.failureMessage = sprintf( ...
            'Only %d fields contain enough signed samples for a joint fit.', ...
            numel(eligibleIndices));
        if ~opts.jointFallbackToIndependent
            error(diagnostics.failureIdentifier, '%s', ...
                diagnostics.failureMessage);
        end
        return
    end

    jointData = repmat(emptyJointFieldData(), 1, numel(eligibleIndices));
    for jointLocalIndex = 1:numel(eligibleIndices)
        globalIndex = eligibleIndices(jointLocalIndex);
        sourceResult = initialResults(globalIndex);
        jointData(jointLocalIndex).globalFieldIndex = globalIndex;
        jointData(jointLocalIndex).t = double(sourceResult.tSeconds(:));
        jointData(jointLocalIndex).timeWeights = jointTimeWeights( ...
            jointData(jointLocalIndex).t, opts);
        diagnostics.timeWeights{globalIndex} = ...
            jointData(jointLocalIndex).timeWeights;
        jointData(jointLocalIndex).y = real(sourceResult.displaySignal(:));
        jointData(jointLocalIndex).isMagnitude = ...
            strcmp(sourceResult.method, 'magnitude_all_zero_crossings');
        if jointData(jointLocalIndex).isMagnitude
            jointData(jointLocalIndex).magnitude = ...
                abs(real(sourceResult.rawSignal(:)));
        else
            jointData(jointLocalIndex).magnitude = [];
        end
        if numel(sourceResult.weights) == numel(sourceResult.tSeconds) && ...
                all(isfinite(sourceResult.weights(:)))
            jointData(jointLocalIndex).robustWeights = min(max( ...
                double(sourceResult.weights(:)), opts.minimumWeight), 1);
        else
            jointData(jointLocalIndex).robustWeights = ...
                ones(size(jointData(jointLocalIndex).t));
        end
        jointData(jointLocalIndex).noiseSigma = sourceResult.noiseSigma;
        jointData(jointLocalIndex).branchIndex = ...
            sourceResult.signChangeIndex;
        jointData(jointLocalIndex).signAmbiguous = ...
            sourceResult.signAmbiguous;
    end

    eligibleRatios = fieldRatios(eligibleIndices);
    initialRates = initialJointRates(initialResults(eligibleIndices), ...
        evolutionFields_mT(eligibleIndices), opts);
    diagnostics.initialR1(eligibleIndices) = initialRates;
    noisePrecision = jointNoisePrecision(jointData, opts);
    diagnostics.noisePrecision(eligibleIndices) = noisePrecision;

    % The independent magnitude fit is deliberately permissive and can
    % initialise a monotonically increasing magnitude trace as all-positive.
    % Before the first joint solve, reconstruct every magnitude branch against
    % a robust estimate of the shared equilibrium scale using the actual
    % inversion constraint S_f(0)<=-alpha_min*M0_D.
    initialSharedM0 = estimateInitialJointM0( ...
        initialResults(eligibleIndices), jointData, eligibleRatios);
    diagnostics.initialSharedM0Detection = initialSharedM0;
    initialBranchWeights = cell(1, numel(jointData));
    for initialBranchIndex = 1:numel(jointData)
        initialBranchWeights{initialBranchIndex} = ...
            jointData(initialBranchIndex).robustWeights;
    end
    [jointData, ~, branchRates] = refineJointMagnitudeBranches( ...
        jointData, eligibleRatios, initialSharedM0, ...
        initialBranchWeights, opts);
    validBranchRates = isfinite(branchRates) & branchRates > 0;
    initialRates(validBranchRates) = branchRates(validBranchRates);

    totalSignPasses = max(1, opts.jointSignIterations);
    jointFit = struct();
    for jointSignPass = 1:totalSignPasses
        jointFit = runJointRobustFit(jointData, eligibleRatios, ...
            initialRates, noisePrecision, opts);
        initialRates = jointFit.rates;
        if jointSignPass >= totalSignPasses
            break
        end
        [jointData, branchesChanged, branchRates] = ...
            refineJointMagnitudeBranches( ...
            jointData, eligibleRatios, jointFit.sharedM0, ...
            jointFit.robustWeights, opts);
        validBranchRates = isfinite(branchRates) & branchRates > 0;
        initialRates(validBranchRates) = branchRates(validBranchRates);
        if ~branchesChanged
            break
        end
    end

    [rateSE, covarianceInfo] = jointRateUncertainty(jointData, ...
        eligibleRatios, jointFit, opts);
    updatedResults = updateJointFieldResults(initialResults, jointData, ...
        eligibleIndices, fieldRatios, jointFit, rateSE, opts);

    diagnostics.used = true;
    diagnostics.finalR1(eligibleIndices) = jointFit.rates;
    diagnostics.sharedM0Detection = jointFit.sharedM0;
    diagnostics.sharedCPer_mT = jointFit.sharedM0 / detectionField_mT;
    diagnostics.inversionAmplitude(eligibleIndices) = ...
        jointFit.inversionAmplitude;
    diagnostics.Q(eligibleIndices) = jointFit.Q;
    diagnostics.alpha(eligibleIndices) = jointFit.alpha;
    diagnostics.weightedSSE = jointFit.weightedSSE;
    diagnostics.globalDfe = covarianceInfo.dfe;
    diagnostics.normalMatrixRcond = covarianceInfo.normalMatrixRcond;
    diagnostics.robustIterations = jointFit.robustIterations;
    diagnostics.rateSweeps = jointFit.rateSweeps;
    diagnostics.optimizerFunctionEvaluations = ...
        jointFit.optimizerFunctionEvaluations;
    diagnostics.converged = jointFit.converged;
    for jointLocalIndex = 1:numel(eligibleIndices)
        diagnostics.magnitudeBranchIndex(eligibleIndices(jointLocalIndex)) = ...
            jointData(jointLocalIndex).branchIndex;
    end
end

function data = emptyJointFieldData()
    data = struct('globalFieldIndex', NaN, 't', [], 'y', [], ...
        'magnitude', [], 'isMagnitude', false, 'robustWeights', [], ...
        'timeWeights', [], 'noiseSigma', NaN, 'branchIndex', NaN, ...
        'signAmbiguous', false);
end

function weights = jointTimeWeights(time, opts)
    time = double(time(:));
    weights = ones(size(time));
    if ~opts.jointTimeWeighting || numel(time) < 2 || ...
            opts.jointLateTimeWeight == 1
        return
    end
    timeMinimum = min(time);
    timeSpan = max(time) - timeMinimum;
    if ~isfinite(timeSpan) || timeSpan <= eps(max(abs(time)))
        return
    end
    normalizedTime = min(max((time - timeMinimum) ./ timeSpan, 0), 1);
    weights = 1 + (opts.jointLateTimeWeight - 1) .* ...
        normalizedTime .^ opts.jointTimeWeightExponent;
    % Keep each field's total leverage unchanged. Only its distribution over
    % evolution time changes, so fields with longer absolute schedules do not
    % dominate fields acquired on shorter schedules.
    weights = weights ./ mean(weights);
end

function rates = initialJointRates(initialResults, eligibleFields_mT, opts)
    rateLower = opts.R1Bounds(1);
    rateUpper = opts.R1Bounds(2);
    rates = nan(1, numel(initialResults));
    reliable = false(size(rates));
    for jointInitialIndex = 1:numel(initialResults)
        value = initialResults(jointInitialIndex).R1;
        if isfinite(value) && value > rateLower && value < rateUpper
            rates(jointInitialIndex) = value;
            reliable(jointInitialIndex) = ...
                ~initialResults(jointInitialIndex).boundHit;
        end
    end
    if any(reliable)
        reliableIndices = find(reliable);
        positiveEligibleFields = eligibleFields_mT(eligibleFields_mT > 0);
        if isempty(positiveEligibleFields)
            fieldFloor = 1;
        else
            fieldFloor = min(positiveEligibleFields) / 10;
        end
        logarithmicFields = log(max(eligibleFields_mT, fieldFloor));
        for missingIndex = find(~reliable)
            [~, nearestPosition] = min(abs( ...
                logarithmicFields(reliableIndices) - ...
                logarithmicFields(missingIndex)));
            rates(missingIndex) = ...
                rates(reliableIndices(nearestPosition));
        end
    else
        rates(:) = sqrt(rateLower * rateUpper);
    end
    rates(~isfinite(rates)) = sqrt(rateLower * rateUpper);
    rates = min(max(rates, rateLower * (1 + 1e-6)), ...
        rateUpper * (1 - 1e-6));
end

function sharedM0 = estimateInitialJointM0( ...
        initialResults, jointData, fieldRatios)
    maximumRatio = max(fieldRatios);
    preferredCandidates = [];
    fallbackCandidates = [];
    for jointM0Index = 1:numel(initialResults)
        ratio = fieldRatios(jointM0Index);
        asymptote = initialResults(jointM0Index).B;
        if ratio <= 0 || ~isfinite(asymptote) || asymptote <= 0
            continue
        end
        candidate = asymptote / ratio;
        if ~isfinite(candidate) || candidate <= 0
            continue
        end
        fallbackCandidates(end + 1) = candidate; %#ok<AGROW>
        if ratio >= 0.25 * maximumRatio && ...
                ~initialResults(jointM0Index).boundHit
            preferredCandidates(end + 1) = candidate; %#ok<AGROW>
        end
    end
    if ~isempty(preferredCandidates)
        candidates = preferredCandidates;
    else
        candidates = fallbackCandidates;
    end
    if ~isempty(candidates)
        centre = median(candidates);
        consistent = candidates >= centre / 10 & candidates <= centre * 10;
        if any(consistent)
            candidates = candidates(consistent);
        end
        sharedM0 = median(candidates);
        return
    end

    % Last-resort scale estimate from the field with the largest equilibrium
    % ratio. It is used only to initialise polarity; the joint fit re-estimates
    % M0_D immediately afterwards.
    [largestRatio, largestIndex] = max(fieldRatios);
    if largestRatio > 0
        observations = abs(real(jointData(largestIndex).y(:)));
        observations = observations(isfinite(observations));
        if ~isempty(observations)
            sharedM0 = max(observations) / largestRatio;
        else
            sharedM0 = 1;
        end
    else
        sharedM0 = 1;
    end
    sharedM0 = max(sharedM0, eps);
end

function precision = jointNoisePrecision(jointData, opts)
    precision = ones(1, numel(jointData));
    if ~opts.jointNoiseWeighting
        return
    end
    sigma = nan(size(precision));
    for jointNoiseIndex = 1:numel(jointData)
        value = jointData(jointNoiseIndex).noiseSigma;
        if isfinite(value) && value > eps
            sigma(jointNoiseIndex) = value;
        end
    end
    validSigma = sigma(isfinite(sigma) & sigma > 0);
    if isempty(validSigma)
        return
    end
    referenceSigma = median(validSigma);
    sigma(~isfinite(sigma) | sigma <= 0) = referenceSigma;
    precision = (referenceSigma ./ sigma).^2;
    precision = min(max(precision, ...
        1 / opts.jointMaximumNoiseWeightRatio), ...
        opts.jointMaximumNoiseWeightRatio);
end

function fit = runJointRobustFit(jointData, fieldRatios, ...
        initialRates, noisePrecision, opts)
    robustWeights = cell(1, numel(jointData));
    for jointWeightIndex = 1:numel(jointData)
        robustWeights{jointWeightIndex} = ...
            jointData(jointWeightIndex).robustWeights(:);
    end
    rates = initialRates;
    previousRates = nan(size(rates));
    previousM0 = NaN;
    totalSweeps = 0;
    totalFunctionEvaluations = 0;
    converged = false;

    for jointRobustIndex = 1:max(1, opts.jointRobustIterations)
        totalWeights = combineJointWeights( ...
            robustWeights, noisePrecision, jointData);
        [rates, linearFit, optimizerInfo] = optimiseJointRates( ...
            jointData, fieldRatios, totalWeights, rates, opts);
        totalSweeps = totalSweeps + optimizerInfo.sweeps;
        totalFunctionEvaluations = totalFunctionEvaluations + ...
            optimizerInfo.functionEvaluations;

        newRobustWeights = cell(size(robustWeights));
        maximumWeightChange = 0;
        for jointResidualIndex = 1:numel(jointData)
            residual = linearFit.residual{jointResidualIndex};
            noiseSigma = robustRealScale(residual);
            if opts.robust
                proposedWeights = huberWeights(abs(residual), ...
                    noiseSigma, opts);
                newRobustWeights{jointResidualIndex} = ...
                    0.5 .* robustWeights{jointResidualIndex} + ...
                    0.5 .* proposedWeights;
            else
                newRobustWeights{jointResidualIndex} = ...
                    ones(size(residual));
            end
            maximumWeightChange = max(maximumWeightChange, max(abs( ...
                newRobustWeights{jointResidualIndex} - ...
                robustWeights{jointResidualIndex})));
        end
        if all(isfinite(previousRates))
            maximumRateChange = max(abs(log(rates ./ previousRates)));
        else
            maximumRateChange = Inf;
        end
        if isfinite(previousM0)
            m0Change = abs(linearFit.sharedM0 - previousM0) / ...
                max(abs(previousM0), eps);
        else
            m0Change = Inf;
        end
        robustWeights = newRobustWeights;
        previousRates = rates;
        previousM0 = linearFit.sharedM0;
        if maximumRateChange < opts.jointTolerance && ...
                m0Change < opts.jointTolerance && ...
                maximumWeightChange < opts.jointTolerance
            converged = true;
            break
        end
    end

    totalWeights = combineJointWeights( ...
        robustWeights, noisePrecision, jointData);
    [rates, linearFit, optimizerInfo] = optimiseJointRates( ...
        jointData, fieldRatios, totalWeights, rates, opts);
    totalSweeps = totalSweeps + optimizerInfo.sweeps;
    totalFunctionEvaluations = totalFunctionEvaluations + ...
        optimizerInfo.functionEvaluations;

    fit = linearFit;
    fit.rates = rates;
    fit.robustWeights = robustWeights;
    fit.totalWeights = totalWeights;
    fit.robustIterations = jointRobustIndex;
    fit.rateSweeps = totalSweeps;
    fit.optimizerFunctionEvaluations = totalFunctionEvaluations;
    fit.converged = converged || optimizerInfo.converged;
end

function totalWeights = combineJointWeights( ...
        robustWeights, noisePrecision, jointData)
    totalWeights = cell(size(robustWeights));
    for jointWeightIndex = 1:numel(robustWeights)
        totalWeights{jointWeightIndex} = ...
            robustWeights{jointWeightIndex} .* ...
            jointData(jointWeightIndex).timeWeights .* ...
            noisePrecision(jointWeightIndex);
    end
end

function [rates, linearFit, optimizerInfo] = optimiseJointRates( ...
        jointData, fieldRatios, totalWeights, initialRates, opts)
    rates = min(max(initialRates(:).', opts.R1Bounds(1)), ...
        opts.R1Bounds(2));
    logLower = log(opts.R1Bounds(1));
    logUpper = log(opts.R1Bounds(2));
    fminOptions = optimset('Display', 'off', 'TolX', 1e-7, ...
        'MaxIter', 100, 'MaxFunEvals', 240);
    functionEvaluations = 0;
    converged = false;

    for jointSweepIndex = 1:max(1, opts.jointRateSweeps)
        previousRates = rates;
        for jointRateIndex = 1:numel(rates)
            objective = @(logRate) jointObjectiveWithRate(logRate, ...
                jointRateIndex, rates, jointData, fieldRatios, ...
                totalWeights, opts);
            [refinedLogRate, refinedValue, ~, output] = fminbnd( ...
                objective, logLower, logUpper, fminOptions);
            functionEvaluations = functionEvaluations + output.funcCount;
            candidateLogRates = [log(rates(jointRateIndex)), ...
                refinedLogRate, logLower, logUpper];
            candidateValues = nan(size(candidateLogRates));
            for jointCandidateIndex = 1:numel(candidateLogRates)
                if jointCandidateIndex == 2
                    candidateValues(jointCandidateIndex) = refinedValue;
                else
                    candidateValues(jointCandidateIndex) = ...
                        objective(candidateLogRates(jointCandidateIndex));
                    functionEvaluations = functionEvaluations + 1;
                end
            end
            [~, bestCandidate] = min(candidateValues);
            rates(jointRateIndex) = ...
                exp(candidateLogRates(bestCandidate));
        end
        if max(abs(log(rates ./ previousRates))) < opts.jointTolerance
            converged = true;
            break
        end
    end
    linearFit = solveJointLinearAtRates(jointData, fieldRatios, ...
        totalWeights, rates, opts);
    optimizerInfo = struct('sweeps', jointSweepIndex, ...
        'functionEvaluations', functionEvaluations, ...
        'converged', converged);
end

function value = jointObjectiveWithRate(logRate, rateIndex, rates, ...
        jointData, fieldRatios, totalWeights, opts)
    trialRates = rates;
    trialRates(rateIndex) = exp(logRate);
    trialFit = solveJointLinearAtRates(jointData, fieldRatios, ...
        totalWeights, trialRates, opts);
    value = trialFit.weightedSSE;
    if ~isfinite(value)
        value = realmax('double');
    end
end

function fit = solveJointLinearAtRates(jointData, fieldRatios, ...
        totalWeights, rates, opts)
    observationCount = sum(arrayfun(@(item) numel(item.t), jointData));
    fieldCount = numel(jointData);
    design = zeros(observationCount, 1 + fieldCount);
    observation = zeros(observationCount, 1);
    squareRootWeight = zeros(observationCount, 1);
    rowStart = 1;
    rowRanges = cell(1, fieldCount);
    for jointDesignIndex = 1:fieldCount
        nObservation = numel(jointData(jointDesignIndex).t);
        rows = rowStart:(rowStart + nObservation - 1);
        rowRanges{jointDesignIndex} = rows;
        exponential = exp(-rates(jointDesignIndex) .* ...
            jointData(jointDesignIndex).t(:));
        % Paper model:
        %   S=M0_D*[r*(1-e)-alpha*e].
        % Write alpha=alpha_min+D/M0_D with D>=0. The linear
        % coefficients are then [M0_D,D_1,...,D_F], all non-negative.
        % This enforces a genuinely inverted initial state and removes the
        % negative-alpha saturation solutions admitted by the v5 Q fit.
        design(rows, 1) = fieldRatios(jointDesignIndex) .* ...
            (1 - exponential) - opts.jointMinimumAlpha .* exponential;
        design(rows, 1 + jointDesignIndex) = -exponential;
        observation(rows) = jointData(jointDesignIndex).y(:);
        squareRootWeight(rows) = sqrt(max( ...
            totalWeights{jointDesignIndex}(:), 0));
        rowStart = rowStart + nObservation;
    end
    weightedDesign = bsxfun(@times, design, squareRootWeight);
    weightedObservation = observation .* squareRootWeight;
    coefficient = nonnegativeCoordinateLeastSquares(weightedDesign, ...
        weightedObservation, opts.jointTolerance);
    predictionVector = design * coefficient;
    residualVector = observation - predictionVector;

    fit = struct();
    fit.sharedM0 = coefficient(1);
    fit.alphaExcessAmplitude = coefficient(2:end).';
    fit.inversionAmplitude = opts.jointMinimumAlpha .* fit.sharedM0 + ...
        fit.alphaExcessAmplitude;
    fit.Q = fieldRatios .* fit.sharedM0 + fit.inversionAmplitude;
    if fit.sharedM0 > eps
        fit.alpha = fit.inversionAmplitude ./ fit.sharedM0;
    else
        fit.alpha = nan(1, fieldCount);
    end
    fit.prediction = cell(1, fieldCount);
    fit.residual = cell(1, fieldCount);
    for jointDesignIndex = 1:fieldCount
        rows = rowRanges{jointDesignIndex};
        fit.prediction{jointDesignIndex} = predictionVector(rows);
        fit.residual{jointDesignIndex} = residualVector(rows);
    end
    fit.weightedSSE = sum(squareRootWeight.^2 .* residualVector.^2);
end

function coefficient = nonnegativeCoordinateLeastSquares( ...
        design, observation, tolerance)
    parameterCount = size(design, 2);
    if isempty(design) || parameterCount == 0
        coefficient = [];
        return
    end
    normalMatrix = design' * design;
    normalRightHandSide = design' * observation;
    conditionEstimate = rcond(normalMatrix);
    if isfinite(conditionEstimate) && conditionEstimate > 1e-12
        unconstrained = real(normalMatrix \ normalRightHandSide);
    else
        unconstrained = real(pinv(normalMatrix) * normalRightHandSide);
    end
    positivityTolerance = max(tolerance, 1e-10) .* ...
        max(norm(unconstrained, inf), 1);
    if all(unconstrained >= -positivityTolerance)
        coefficient = max(unconstrained, 0);
        return
    end
    coefficient = max(unconstrained, 0);
    columnEnergy = sum(design.^2, 1);
    residual = observation - design * coefficient;
    maximumIterations = 200;
    for nnlsIteration = 1:maximumIterations
        previousCoefficient = coefficient;
        for nnlsIndex = 1:parameterCount
            if columnEnergy(nnlsIndex) <= eps
                coefficient(nnlsIndex) = 0;
                continue
            end
            oldCoefficient = coefficient(nnlsIndex);
            partialResidual = residual + ...
                design(:, nnlsIndex) .* oldCoefficient;
            newCoefficient = max(0, ...
                design(:, nnlsIndex)' * partialResidual / ...
                columnEnergy(nnlsIndex));
            coefficient(nnlsIndex) = newCoefficient;
            residual = partialResidual - ...
                design(:, nnlsIndex) .* newCoefficient;
        end
        relativeChange = norm(coefficient - previousCoefficient) / ...
            max(norm(previousCoefficient), 1);
        if relativeChange < max(tolerance * 0.1, 1e-8)
            break
        end
    end
end

function [jointData, changed, selectedRates] = refineJointMagnitudeBranches( ...
        jointData, fieldRatios, sharedM0, robustWeights, opts)
    changed = false;
    selectedRates = nan(1, numel(jointData));
    for jointBranchFieldIndex = 1:numel(jointData)
        if ~jointData(jointBranchFieldIndex).isMagnitude
            continue
        end
        magnitude = jointData(jointBranchFieldIndex).magnitude(:);
        time = jointData(jointBranchFieldIndex).t(:);
        intercept = sharedM0 * fieldRatios(jointBranchFieldIndex);
        weights = robustWeights{jointBranchFieldIndex}(:) .* ...
            jointData(jointBranchFieldIndex).timeWeights(:);
        sampleCount = numel(time);
        candidateScore = inf(1, sampleCount + 1);
        candidateRate = nan(1, sampleCount + 1);
        candidateQ = nan(1, sampleCount + 1);
        for branchIndex = 0:sampleCount
            signs = ones(sampleCount, 1);
            if branchIndex > 0
                signs(1:branchIndex) = -1;
            end
            signedSignal = signs .* magnitude;
            [candidateRate(branchIndex + 1), ...
                candidateQ(branchIndex + 1), ...
                candidateScore(branchIndex + 1)] = ...
                fitFixedInterceptBranch(time, signedSignal, weights, ...
                intercept, sharedM0, opts);
        end
        [sortedScores, order] = sort(candidateScore, 'ascend');
        if isempty(order) || ~isfinite(sortedScores(1))
            continue
        end
        bestBranch = order(1) - 1;
        signs = ones(sampleCount, 1);
        if bestBranch > 0
            signs(1:bestBranch) = -1;
        end
        newSignal = signs .* magnitude;
        if numel(newSignal) ~= numel(jointData(jointBranchFieldIndex).y) || ...
                any(newSignal ~= jointData(jointBranchFieldIndex).y)
            changed = true;
            jointData(jointBranchFieldIndex).y = newSignal;
        end
        jointData(jointBranchFieldIndex).branchIndex = bestBranch;
        selectedRates(jointBranchFieldIndex) = candidateRate(order(1));
        if numel(sortedScores) > 1 && isfinite(sortedScores(2))
            signalEnergy = max(sum(magnitude.^2), eps);
            closeScore = sortedScores(2) <= opts.signAmbiguityRatio * ...
                max(sortedScores(1), eps * signalEnergy);
            differentRate = abs(log(candidateRate(order(2)) / ...
                candidateRate(order(1)))) > ...
                opts.signAmbiguityLogRateDifference;
            jointData(jointBranchFieldIndex).signAmbiguous = ...
                closeScore && differentRate;
        else
            jointData(jointBranchFieldIndex).signAmbiguous = false;
        end
        minimumQ = intercept + opts.jointMinimumAlpha * sharedM0;
        if ~isfinite(candidateQ(order(1))) || ...
                candidateQ(order(1)) < minimumQ
            jointData(jointBranchFieldIndex).signAmbiguous = true;
        end
    end
end

function [rate, Q, score] = fitFixedInterceptBranch( ...
        time, signal, weights, intercept, sharedM0, opts)
    logLower = log(opts.R1Bounds(1));
    logUpper = log(opts.R1Bounds(2));
    minimumQ = intercept + opts.jointMinimumAlpha * sharedM0;
    objective = @(logRate) fixedInterceptProfile(logRate, time, ...
        signal, weights, intercept, minimumQ);
    fminOptions = optimset('Display', 'off', 'TolX', 1e-8, ...
        'MaxIter', 100, 'MaxFunEvals', 220);
    refinedLogRate = fminbnd(objective, logLower, logUpper, fminOptions);
    candidateLogRate = [refinedLogRate, logLower, logUpper];
    candidateScore = arrayfun(objective, candidateLogRate);
    [score, bestIndex] = min(candidateScore);
    rate = exp(candidateLogRate(bestIndex));
    [Q, score] = fixedInterceptLinearFit(rate, time, signal, ...
        weights, intercept, minimumQ);
    if Q < minimumQ
        score = Inf;
        return
    end
    tolerance = 0.01 * (logUpper - logLower);
    if candidateLogRate(bestIndex) <= logLower + tolerance || ...
            candidateLogRate(bestIndex) >= logUpper - tolerance
        score = score + opts.boundCandidatePenaltyFraction * ...
            max(sum(signal.^2), eps);
    end
end

function score = fixedInterceptProfile(logRate, time, signal, ...
        weights, intercept, minimumQ)
    [~, score] = fixedInterceptLinearFit(exp(logRate), time, ...
        signal, weights, intercept, minimumQ);
end

function [Q, score] = fixedInterceptLinearFit(rate, time, signal, ...
        weights, intercept, minimumQ)
    exponential = exp(-rate .* time(:));
    weights = max(weights(:), 0);
    denominator = sum(weights .* exponential.^2);
    if denominator <= eps
        Q = minimumQ;
    else
        Q = max(sum(weights .* exponential .* ...
            (intercept - signal(:))) / denominator, minimumQ);
    end
    residual = signal(:) - (intercept - Q .* exponential);
    score = sum(weights .* residual.^2);
end

function [rateSE, info] = jointRateUncertainty(jointData, ...
        fieldRatios, jointFit, opts)
    fieldCount = numel(jointData);
    observationCount = sum(arrayfun(@(item) numel(item.t), jointData));
    parameterCount = 1 + 2 * fieldCount;
    jacobian = zeros(observationCount, parameterCount);
    residual = zeros(observationCount, 1);
    squareRootWeight = zeros(observationCount, 1);
    rowStart = 1;
    for jointJacobianIndex = 1:fieldCount
        time = jointData(jointJacobianIndex).t(:);
        nObservation = numel(time);
        rows = rowStart:(rowStart + nObservation - 1);
        exponential = exp(-jointFit.rates(jointJacobianIndex) .* time);
        jacobian(rows, 1) = fieldRatios(jointJacobianIndex);
        jacobian(rows, 1 + jointJacobianIndex) = -exponential;
        jacobian(rows, 1 + fieldCount + jointJacobianIndex) = ...
            jointFit.Q(jointJacobianIndex) .* time .* exponential;
        residual(rows) = jointFit.residual{jointJacobianIndex};
        squareRootWeight(rows) = sqrt(max( ...
            jointFit.totalWeights{jointJacobianIndex}(:), 0));
        rowStart = rowStart + nObservation;
    end
    weightedJacobian = bsxfun(@times, jacobian, squareRootWeight);
    normalMatrix = weightedJacobian' * weightedJacobian;
    dfe = max(observationCount - parameterCount, 1);
    varianceEstimate = max(sum((squareRootWeight .* residual).^2) / dfe, eps);
    covariance = pinv(normalMatrix) .* varianceEstimate;
    rateVariance = diag(covariance( ...
        (2 + fieldCount):(1 + 2 * fieldCount), ...
        (2 + fieldCount):(1 + 2 * fieldCount)));
    rateSE = sqrt(max(real(rateVariance(:).'), 0));
    rateSE(~isfinite(rateSE)) = NaN;
    info = struct('dfe', dfe, 'normalMatrixRcond', rcond(normalMatrix));
end

function updatedResults = updateJointFieldResults(initialResults, ...
        jointData, eligibleIndices, allFieldRatios, jointFit, rateSE, opts)
    updatedResults = initialResults;
    zCritical = sqrt(2) * erfinv(opts.confidenceLevel);
    logSpan = log(opts.R1Bounds(2)) - log(opts.R1Bounds(1));
    for jointUpdateIndex = 1:numel(eligibleIndices)
        globalIndex = eligibleIndices(jointUpdateIndex);
        initialResult = initialResults(globalIndex);
        result = initialResult;
        time = jointData(jointUpdateIndex).t(:);
        signal = jointData(jointUpdateIndex).y(:);
        prediction = jointFit.prediction{jointUpdateIndex};
        residual = jointFit.residual{jointUpdateIndex};
        rate = jointFit.rates(jointUpdateIndex);
        Q = jointFit.Q(jointUpdateIndex);
        intercept = jointFit.sharedM0 * allFieldRatios(globalIndex);
        sigma = robustRealScale(residual);
        sse = sum(residual.^2);
        sst = sum((signal - mean(signal)).^2);
        dfe = max(numel(time) - 2, 1);
        rmse = sqrt(sse / dfe);
        scaleAmplitude = max([Q, ...
            max(abs(signal - median(signal))), eps]);

        result.ok = isfinite(rate) && rate > 0 && ...
            isfinite(Q) && isfinite(intercept);
        result.method = 'joint_multifield_inversion_constrained';
        result.tSeconds = time;
        result.displaySignal = signal;
        result.displayFit = prediction;
        result.predictionComplex = [];
        result.residual = residual;
        result.weights = jointFit.robustWeights{jointUpdateIndex};
        result.outlierIndices = find( ...
            result.weights < opts.outlierWeightThreshold).';
        result.B = intercept;
        result.C = Q;
        result.S0 = intercept - Q;
        result.Sinf = intercept;
        result.R1 = rate;
        result.T1ms = 1000 / rate;
        if isfinite(rateSE(jointUpdateIndex))
            rawRateCI = rate + [-1 1] .* ...
                zCritical .* rateSE(jointUpdateIndex);
            result.ciOpen = [rawRateCI(1) <= opts.R1Bounds(1), ...
                rawRateCI(2) >= opts.R1Bounds(2)];
            result.R1CI95 = [max(rawRateCI(1), opts.R1Bounds(1)), ...
                min(rawRateCI(2), opts.R1Bounds(2))];
            if result.R1CI95(1) <= 0
                result.R1CI95(1) = opts.R1Bounds(1);
            end
            result.T1CI95 = rateCiToT1Ci(result.R1CI95);
            result.R1SE = rateSE(jointUpdateIndex);
            result.T1SE = 1000 .* rateSE(jointUpdateIndex) ./ rate.^2;
            if any(result.ciOpen)
                % A symmetric covariance error is misleading when the
                % interval reaches a rate bound and can dominate the entire
                % dispersion plot by several orders of magnitude.
                result.R1SE = NaN;
                result.T1SE = NaN;
            end
        else
            result.R1CI95 = [NaN NaN];
            result.T1CI95 = [NaN NaN];
            result.R1SE = NaN;
            result.T1SE = NaN;
            result.ciOpen = [false false];
        end
        result.betaComplex = [intercept, -Q];
        result.orthogonalOffset = NaN;
        result.orthogonalNRMSE = NaN;
        result.noiseSigma = sigma;
        result.dynamicSNR = Q / max(sigma, eps);
        result.sse = sse;
        result.rsquare = safeRsquare(sse, sst);
        result.rmse = rmse;
        result.nrmse = rmse / scaleAmplitude;
        result.dfe = dfe;
        result.boundHit = abs(log(rate / opts.R1Bounds(1))) <= ...
            0.01 * logSpan || abs(log(opts.R1Bounds(2) / rate)) <= ...
            0.01 * logSpan;
        if result.boundHit
            result.R1SE = NaN;
            result.T1SE = NaN;
        end
        result.zeroCrossingSeconds = ...
            zeroCrossingTime(intercept, Q, rate);
        result.signChangeIndex = sum(signal < 0);
        result.reliableSignTransitions = countReliableSignTransitions( ...
            signal, sigma, Q);
        result.signAmbiguous = jointData(jointUpdateIndex).signAmbiguous;
        result.jointSharedM0 = jointFit.sharedM0;
        result.jointFieldRatio = allFieldRatios(globalIndex);
        result.jointAlpha = jointFit.alpha(jointUpdateIndex);
        result.jointTimeWeights = ...
            jointData(jointUpdateIndex).timeWeights;
        result.jointTotalWeights = ...
            jointFit.totalWeights{jointUpdateIndex};
        result.independentFitSummary = fitSummary(initialResult);
        if isempty(initialResult.selectionReason)
            result.selectionReason = ...
                'joint all-field inversion-constrained refinement';
        else
            result.selectionReason = [initialResult.selectionReason, ...
                '; joint all-field inversion-constrained refinement'];
        end
        result.qc = 'not fitted';
        result.qcReasons = {};
        result = applyQualityControl(result, opts);
        if ~jointFit.converged && ~strncmp(result.qc, 'failed', 6)
            result.qcReasons{end + 1} = ...
                'joint optimizer did not meet the requested tolerance';
            result.qc = ['review: ' strjoin(result.qcReasons, '; ')];
        end
        updatedResults(globalIndex) = result;
    end
end


% ========================================================================
% Variable projection, robust weights and profile intervals
% ========================================================================

function [R1, beta, sse, boundHit] = optimiseRate(t, y, weights, opts)
    logLower = log(opts.R1Bounds(1));
    logUpper = log(opts.R1Bounds(2));
    logGrid = linspace(logLower, logUpper, opts.rateGridPoints);
    gridSse = nan(size(logGrid));
    for k = 1:numel(logGrid)
        gridSse(k) = profileSseAtLogRate(logGrid(k), t, y, weights);
    end
    [~, bestIndex] = min(gridSse);

    leftIndex = max(1, bestIndex - 1);
    rightIndex = min(numel(logGrid), bestIndex + 1);
    localLower = logGrid(leftIndex);
    localUpper = logGrid(rightIndex);
    candidates = [logGrid(bestIndex), logLower, logUpper];
    if localUpper > localLower
        fminOptions = optimset('Display', 'off', 'TolX', 1e-8, ...
            'MaxIter', 100, 'MaxFunEvals', 200);
        refined = fminbnd(@(lr) profileSseAtLogRate(lr, t, y, weights), ...
            localLower, localUpper, fminOptions);
        candidates(end + 1) = refined; %#ok<AGROW>
    end

    candidateSse = nan(size(candidates));
    for k = 1:numel(candidates)
        candidateSse(k) = profileSseAtLogRate(candidates(k), t, y, weights);
    end
    [sse, selected] = min(candidateSse);
    selectedLogRate = candidates(selected);
    R1 = exp(selectedLogRate);
    [beta, sse] = solveLinearAtRate(t, y, weights, R1);

    tolerance = 0.01 * (logUpper - logLower);
    boundHit = selectedLogRate <= logLower + tolerance || ...
        selectedLogRate >= logUpper - tolerance;
end

function sse = profileSseAtLogRate(logRate, t, y, weights)
    R1 = exp(logRate);
    [~, sse] = solveLinearAtRate(t, y, weights, R1);
    if ~isfinite(sse)
        sse = realmax('double');
    end
end

function [beta, sse] = solveLinearAtRate(t, y, weights, R1)
    exponential = exp(-R1 .* t(:));
    X = [ones(size(exponential)), exponential];
    squareRootWeight = sqrt(max(weights(:), 0));
    Xw = bsxfun(@times, X, squareRootWeight);
    yw = y(:) .* squareRootWeight;
    beta = Xw \ yw;
    residual = y(:) - X * beta;
    sse = sum(weights(:) .* abs(residual).^2);
end

function prediction = affinePrediction(t, R1, beta)
    prediction = beta(1) + beta(2) .* exp(-R1 .* t(:));
end

function weights = huberWeights(residualMagnitude, noiseSigma, opts)
    residualMagnitude = abs(residualMagnitude(:));
    if ~isfinite(noiseSigma) || noiseSigma <= eps
        weights = ones(size(residualMagnitude));
        return
    end
    cutoff = opts.robustTune * noiseSigma;
    weights = ones(size(residualMagnitude));
    large = residualMagnitude > cutoff;
    weights(large) = cutoff ./ residualMagnitude(large);
    weights = max(weights, opts.minimumWeight);
end

function scale = robustRealScale(residual)
    residual = real(residual(:));
    residual = residual(isfinite(residual));
    if isempty(residual)
        scale = NaN;
        return
    end
    centre = median(residual);
    scale = 1.4826 * median(abs(residual - centre));
    if ~isfinite(scale) || scale <= eps
        scale = sqrt(mean((residual - mean(residual)).^2));
    end
    if ~isfinite(scale) || scale <= eps
        scale = eps;
    end
end

function scale = robustComplexScale(residual)
    residual = residual(:);
    good = isfinite(real(residual)) & isfinite(imag(residual));
    residual = residual(good);
    if isempty(residual)
        scale = NaN;
        return
    end
    scaleReal = robustRealScale(real(residual));
    scaleImag = robustRealScale(imag(residual));
    scale = sqrt(scaleReal.^2 + scaleImag.^2);
    if ~isfinite(scale) || scale <= eps
        scale = sqrt(mean(abs(residual - mean(residual)).^2));
    end
    if ~isfinite(scale) || scale <= eps
        scale = eps;
    end
end

function [rateCI, openBound] = profileRateCI(t, y, weights, R1hat, ...
        realObservationCount, parameterCount, opts)
    logLower = log(opts.R1Bounds(1));
    logUpper = log(opts.R1Bounds(2));
    logGrid = linspace(logLower, logUpper, opts.profileGridPoints);
    logGrid = unique(sort([logGrid, log(R1hat)]));
    profile = nan(size(logGrid));
    for k = 1:numel(logGrid)
        profile(k) = profileSseAtLogRate(logGrid(k), t, y, weights);
    end
    [minimumSse, minimumIndex] = min(profile);
    dfe = max(realObservationCount - parameterCount, 1);
    varianceEstimate = max(minimumSse / dfe, eps);
    chiSquareOneDf = 2 * erfinv(opts.confidenceLevel).^2;
    threshold = minimumSse + chiSquareOneDf * varianceEstimate;
    inside = profile <= threshold;

    leftIndex = minimumIndex;
    while leftIndex > 1 && inside(leftIndex - 1)
        leftIndex = leftIndex - 1;
    end
    rightIndex = minimumIndex;
    while rightIndex < numel(logGrid) && inside(rightIndex + 1)
        rightIndex = rightIndex + 1;
    end

    leftOpen = leftIndex == 1;
    rightOpen = rightIndex == numel(logGrid);
    if leftOpen
        leftLog = logGrid(1);
    else
        leftLog = interpolateThreshold(logGrid(leftIndex - 1), ...
            profile(leftIndex - 1), logGrid(leftIndex), ...
            profile(leftIndex), threshold);
    end
    if rightOpen
        rightLog = logGrid(end);
    else
        rightLog = interpolateThreshold(logGrid(rightIndex), ...
            profile(rightIndex), logGrid(rightIndex + 1), ...
            profile(rightIndex + 1), threshold);
    end
    rateCI = sort(exp([leftLog, rightLog]));
    openBound = [leftOpen, rightOpen];
end

function x = interpolateThreshold(x1, y1, x2, y2, threshold)
    if ~isfinite(y1) || ~isfinite(y2) || y2 == y1
        x = 0.5 * (x1 + x2);
        return
    end
    fraction = (threshold - y1) / (y2 - y1);
    fraction = min(max(fraction, 0), 1);
    x = x1 + fraction * (x2 - x1);
end


% ========================================================================
% Diagnostics, storage and plotting
% ========================================================================

function r = applyQualityControl(r, opts)
    reasons = {};
    if ~r.ok || ~isfinite(r.R1) || r.R1 <= 0
        r.qc = 'failed: non-finite rate';
        r.qcReasons = {r.qc};
        return
    end
    if r.boundHit
        reasons{end + 1} = 'rate at search bound'; %#ok<AGROW>
    end
    if any(r.ciOpen)
        reasons{end + 1} = 'profile interval reaches a search bound'; %#ok<AGROW>
    end
    if all(isfinite(r.T1CI95)) && min(r.T1CI95) > 0 && ...
            max(r.T1CI95) / min(r.T1CI95) > opts.maximumT1CIRatio
        reasons{end + 1} = 'wide T1 profile interval'; %#ok<AGROW>
    end
    if ~isfinite(r.dynamicSNR) || r.dynamicSNR < opts.minimumDynamicSNR
        reasons{end + 1} = 'weak fitted recovery amplitude'; %#ok<AGROW>
    end
    if strncmp(r.method, 'joint_multifield', 16) && ...
            (~isfinite(r.C) || r.C <= 0)
        reasons{end + 1} = ...
            'joint fit reached a zero recovery amplitude'; %#ok<AGROW>
    end
    if strcmp(r.method, 'joint_multifield_inversion_constrained') && ...
            isfinite(r.jointAlpha) && ...
            r.jointAlpha <= opts.jointMinimumAlpha * (1 + 1e-3)
        reasons{end + 1} = ...
            'inversion efficiency reached its lower constraint'; %#ok<AGROW>
    end
    if ~isfinite(r.nrmse) || r.nrmse > opts.maximumNRMSE
        reasons{end + 1} = 'large model residuals'; %#ok<AGROW>
    end
    if any(strcmp(r.method, {'complex_affine_variable_projection', ...
            'complex_phase_drift_ir'})) && ...
            isfinite(r.orthogonalNRMSE) && ...
            r.orthogonalNRMSE > opts.maximumOrthogonalNRMSE
        reasons{end + 1} = 'complex trajectory is not approximately one-dimensional'; %#ok<AGROW>
    end
    if any(strcmp(r.method, {'complex_affine_variable_projection', ...
            'complex_phase_drift_ir'})) && ...
            isfinite(r.reliableSignTransitions) && ...
            r.reliableSignTransitions > opts.maximumReliableSignTransitions
        reasons{end + 1} = 'complex projection has more than one reliable sign transition'; %#ok<AGROW>
    end
    if strcmp(r.method, 'complex_phase_drift_ir') && ...
            (~isfinite(r.phaseResidualDegrees) || ...
             r.phaseResidualDegrees > opts.maximumPhaseResidualDegrees)
        reasons{end + 1} = 'complex phase drift is not sufficiently coherent'; %#ok<AGROW>
    end
    if strcmp(r.method, 'magnitude_rician_ir') && ...
            strcmp(r.noiseSource, 'recovery minimum fallback')
        reasons{end + 1} = ...
            'magnitude noise was estimated from the recovery minimum'; %#ok<AGROW>
    end
    if strcmp(r.method, 'magnitude_rician_ir') && ...
            strcmp(r.magnitudeNoiseModel, 'noncentral_chi_approx')
        reasons{end + 1} = ...
            'multi-coil magnitude noise model is approximate'; %#ok<AGROW>
    end
    if numel(r.outlierIndices) > max(1, floor(0.25 * numel(r.tSeconds)))
        reasons{end + 1} = 'several low-weight time points'; %#ok<AGROW>
    end
    if r.signAmbiguous
        reasons{end + 1} = 'magnitude polarity is ambiguous'; %#ok<AGROW>
    end

    T1seconds = r.T1ms / 1000;
    finiteTimes = r.tSeconds(isfinite(r.tSeconds) & r.tSeconds >= 0);
    if isempty(finiteTimes)
        minimumTime = 0;
    else
        minimumTime = min(finiteTimes);
    end
    maximumTime = max(r.tSeconds);
    r.timeCoverageT1 = [minimumTime / T1seconds, maximumTime / T1seconds];
    if maximumTime / T1seconds < opts.minimumLateCoverageT1
        reasons{end + 1} = 'no late-time leverage on the asymptote'; %#ok<AGROW>
    end
    if minimumTime > 0 && minimumTime / T1seconds > opts.maximumEarlyCoverageT1
        reasons{end + 1} = 'no early-time leverage on the initial state'; %#ok<AGROW>
    end

    r.qcReasons = reasons;
    if isempty(reasons)
        r.qc = 'good';
    else
        r.qc = ['review: ' strjoin(reasons, '; ')];
    end
end

function results = assembleResults(fieldResult, fields_mT, slice, signalMode, opts)
    resultFieldCount = numel(fieldResult);
    results = struct();
    results.modelEquation = 'S(t) = S_inf + (S_0-S_inf)*exp(-R1*t)';
    results.observationEquation = ...
        'complex: phase-drift-corrected S(t); magnitude: E[|S(t)+noise|]';
    results.modelName = 'phase_aware_constrained_field_cycling_IR';
    results.fitMode = opts.fitMode;
    results.signalMode = signalMode;
    results.slice = slice;
    results.fields_mT = fields_mT(:).';
    results.options = opts;
    results.field = fieldResult;
    results.R1 = nan(1, resultFieldCount);
    results.R1SE = nan(1, resultFieldCount);
    results.T1ms = nan(1, resultFieldCount);
    results.T1SE = nan(1, resultFieldCount);
    results.R1CI95 = nan(2, resultFieldCount);
    results.T1CI95 = nan(2, resultFieldCount);
    results.qc = cell(1, resultFieldCount);
    for n = 1:resultFieldCount
        results.R1(n) = fieldResult(n).R1;
        results.R1SE(n) = fieldResult(n).R1SE;
        results.T1ms(n) = fieldResult(n).T1ms;
        results.T1SE(n) = fieldResult(n).T1SE;
        results.R1CI95(:, n) = fieldResult(n).R1CI95(:);
        results.T1CI95(:, n) = fieldResult(n).T1CI95(:);
        results.qc{n} = fieldResult(n).qc;
    end
end

function obj = storeResults(obj, results)
    resultFieldCount = numel(results.field);
    fitParams = cell(1, resultFieldCount);
    fitSse = nan(1, resultFieldCount);
    fitRsq = nan(1, resultFieldCount);
    fitRmse = nan(1, resultFieldCount);
    fitDfe = nan(1, resultFieldCount);
    fitWeights = cell(1, resultFieldCount);
    autoPolarityFlips = cell(1, resultFieldCount);
    phase = nan(1, resultFieldCount);
    phaseTrend = cell(1, resultFieldCount);
    for n = 1:resultFieldCount
        storedFieldResult = results.field(n);
        fitParams{n} = [storedFieldResult.B, storedFieldResult.C, storedFieldResult.R1];
        fitSse(n) = storedFieldResult.sse;
        fitRsq(n) = storedFieldResult.rsquare;
        fitRmse(n) = storedFieldResult.rmse;
        fitDfe(n) = storedFieldResult.dfe;
        fitWeights{n} = storedFieldResult.weights;
        phase(n) = storedFieldResult.phaseRad;
        phaseTrend{n} = storedFieldResult.phaseTrendRad;
        magnitudeDerived = any(strcmp(storedFieldResult.method, ...
            {'magnitude_all_zero_crossings', 'magnitude_rician_ir'})) || ...
            (strcmp(storedFieldResult.method, ...
                'joint_multifield_inversion_constrained') && ...
             isstruct(storedFieldResult.independentFitSummary) && ...
             isfield(storedFieldResult.independentFitSummary, 'method') && ...
             any(strcmp(storedFieldResult.independentFitSummary.method, ...
                {'magnitude_all_zero_crossings', 'magnitude_rician_ir'})));
        if magnitudeDerived && ...
                numel(storedFieldResult.displaySignal) == ...
                numel(storedFieldResult.tSeconds)
            autoPolarityFlips{n} = find(storedFieldResult.displaySignal < 0).';
        else
            autoPolarityFlips{n} = [];
        end
    end

    obj = safeSet(obj, 'dispersioncurve', results.R1);
    obj = safeSet(obj, 'T1FitResults', results);
    obj = safeSet(obj, 'R1dispersion', results.R1);
    obj = safeSet(obj, 'R1error', results.R1SE);
    obj = safeSet(obj, 'T1dispersion', results.T1ms);
    obj = safeSet(obj, 'T1error', results.T1SE);
    obj = safeSet(obj, 'R1ci95', results.R1CI95);
    obj = safeSet(obj, 'T1ci95', results.T1CI95);
    obj = safeSet(obj, 'fit_params', fitParams);
    obj = safeSet(obj, 'fit_sse', fitSse);
    obj = safeSet(obj, 'fit_rsq', fitRsq);
    obj = safeSet(obj, 'fit_rmse', fitRmse);
    obj = safeSet(obj, 'fit_dfe', fitDfe);
    obj = safeSet(obj, 'fit_qc', results.qc);
    obj = safeSet(obj, 'fit_weights', fitWeights);
    obj = safeSet(obj, 'fit_model_name', results.modelName);
    obj = safeSet(obj, 'phase_reference_rad', phase);
    obj = safeSet(obj, 'phase_drift_rad', phaseTrend);
    obj = safeSet(obj, 'fit_diagnostics', results.field);
    % Existing manual polarity and point-edit members are intentionally not
    % reset here. They are user decisions and must survive refitting.
    obj = safeSet(obj, 'auto_polarity_flips', autoPolarityFlips);
end

function printFitSummary(results)
    fprintf('\nRobust field-cycling T1 fit: %s\n', results.signalMode);
    if isfield(results, 'implementation')
        fprintf('Implementation: %s\n', results.implementation);
    end
    if isfield(results, 'inputDiagnostics')
        d = results.inputDiagnostics;
        fprintf(['Input geometry: %d fields | timepoints %s | magimage %s | ', ...
                 'compleximage %s\n'], d.fieldCount, mat2str(d.timepointsSize), ...
            mat2str(d.magimageSize), mat2str(d.compleximageSize));
    end
    fprintf('Model: %s\n', results.modelEquation);
    if isfield(results, 'joint') && results.joint.used
        fprintf(['Joint inversion fit: shared M0_D %.6g | B_D %.6g mT | ', ...
                 'alpha >= %.4g | late/early weight %.4g | ', ...
                 'weighted SSE %.6g\n'], ...
            results.joint.sharedM0Detection, ...
            results.joint.detectionField_mT, results.joint.minimumAlpha, ...
            results.joint.lateTimeWeight, ...
            results.joint.weightedSSE);
    elseif isfield(results, 'joint') && results.joint.attempted && ...
            ~isempty(results.joint.failureMessage)
        fprintf('Joint fit not used: %s\n', results.joint.failureMessage);
    end
    fprintf('%6s %10s %12s %12s %-34s %s\n', ...
        'Index', 'Field_mT', 'R1_s-1', 'T1_ms', 'Method', 'QC');
    for n = 1:numel(results.field)
        summaryFieldResult = results.field(n);
        fprintf('%6d %10.5g %12.5g %12.5g %-34s %s\n', ...
            n, summaryFieldResult.field_mT, summaryFieldResult.R1, ...
            summaryFieldResult.T1ms, summaryFieldResult.method, ...
            summaryFieldResult.qc);
    end
end

function makeDispersionPlot(results, obj)
    valid = isfinite(results.R1) & results.R1 > 0 & isfinite(results.fields_mT);
    if ~any(valid)
        return
    end
    plotR1 = memberExists(obj, 'R1T1') && isscalar(obj.R1T1) && obj.R1T1 == 1;
    figure('Name', 'Robust field-cycling T1 dispersion');
    if plotR1
        x = results.fields_mT(valid) .* 42.57747892e6 ./ 1e6 ./ 1000;
        y = results.R1(valid);
        intervals = results.R1CI95(:, valid);
        lowerError = max(y - intervals(1, :), 0);
        upperError = max(intervals(2, :) - y, 0);
        lowerError(~isfinite(lowerError)) = 0;
        upperError(~isfinite(upperError)) = 0;
        errorbar(x, y, lowerError, upperError, 'o-', 'LineWidth', 1.2);
        xlabel('Evolution field (MHz)');
        ylabel('R_1 (s^{-1})');
    else
        x = results.fields_mT(valid) ./ 1000;
        y = results.T1ms(valid);
        intervals = results.T1CI95(:, valid);
        lowerError = max(y - intervals(1, :), 0);
        upperError = max(intervals(2, :) - y, 0);
        lowerError(~isfinite(lowerError)) = 0;
        upperError(~isfinite(upperError)) = 0;
        errorbar(x, y, lowerError, upperError, 'o-', 'LineWidth', 1.2);
        xlabel('Evolution field (T)');
        ylabel('T_1 (ms)');
    end
    hold on;
    validIndices = find(valid);
    review = ~arrayfun(@(item) strncmp(item.qc, 'good', 4), ...
        results.field(valid));
    if any(review)
        plot(x(review), y(review), 'rx', 'MarkerSize', 9, ...
            'LineWidth', 1.5, 'DisplayName', 'review');
    end
    for labelIndex = 1:numel(validIndices)
        text(x(labelIndex), y(labelIndex), ...
            sprintf('  %d', validIndices(labelIndex)), 'FontSize', 8);
    end
    hold off;
    if all(x > 0)
        set(gca, 'XScale', 'log');
    end
    grid on;
    if isfield(results, 'joint') && results.joint.used
        if results.joint.timeWeightingEnabled
            title(sprintf(['Joint multi-field inversion fit; ', ...
                'late/early weight %.3g'], results.joint.lateTimeWeight));
        else
            title('Joint multi-field inversion-constrained recovery fit');
        end
    else
        title('Independent phase-aware constrained IR fits');
    end
end

function makeQaDashboard(results)
    panelIndex = find(arrayfun(@(r) ~isempty(r.tSeconds), results.field));
    if isempty(panelIndex)
        return
    end
    nPanels = numel(panelIndex);
    nColumns = max(2, ceil(sqrt(nPanels)));
    fitRows = ceil(nPanels / nColumns);
    overviewRows = 2;
    nRows = overviewRows + fitRows;
    figure('Name', 'Robust field-cycling T1 fit QA');
    layout = tiledlayout(nRows, nColumns, 'TileSpacing', 'compact', 'Padding', 'compact');

    overview = nexttile(layout, 1, [overviewRows nColumns]);
    valid = isfinite(results.R1) & results.R1 > 0 & isfinite(results.fields_mT);
    hold(overview, 'on');
    if any(valid)
        yValues = results.T1ms(valid);
        intervals = results.T1CI95(:, valid);
        lowerError = max(yValues - intervals(1, :), 0);
        upperError = max(intervals(2, :) - yValues, 0);
        lowerError(~isfinite(lowerError)) = 0;
        upperError(~isfinite(upperError)) = 0;
        errorbar(overview, results.fields_mT(valid), results.T1ms(valid), ...
            lowerError, upperError, 'o-', 'LineWidth', 1.3, ...
            'DisplayName', 'all fitted fields');
        for n = find(valid)
            text(overview, results.fields_mT(n), results.T1ms(n), ...
                sprintf('  %d', n), 'FontSize', 8);
        end
        reviewIndices = find(valid & ~arrayfun(@(item) ...
            strncmp(item.qc, 'good', 4), results.field));
        if ~isempty(reviewIndices)
            plot(overview, results.fields_mT(reviewIndices), ...
                results.T1ms(reviewIndices), 'rx', 'MarkerSize', 9, ...
                'LineWidth', 1.4, 'DisplayName', 'review');
        end
    end
    if all(results.fields_mT(valid) > 0)
        set(overview, 'XScale', 'log');
    end
    grid(overview, 'on');
    xlabel(overview, 'Evolution field (mT)');
    ylabel(overview, 'T_1 (ms)');
    if isfield(results, 'joint') && results.joint.used
        overviewTitle = sprintf( ...
            ['Joint all-field fit: %d/%d finite T1 | ', ...
             'late/early weight %.3g'], ...
            sum(valid), numel(results.field), results.joint.lateTimeWeight);
    else
        overviewTitle = sprintf( ...
            'Independent constrained IR fits: %d/%d fields returned finite T1', ...
            sum(valid), numel(results.field));
    end
    title(overview, overviewTitle);
    hold(overview, 'off');

    firstFitTile = overviewRows * nColumns + 1;
    for k = 1:nPanels
        dashboardFieldResult = results.field(panelIndex(k));
        ax = nexttile(layout, firstFitTile + k - 1);
        hold(ax, 'on');
        displayed = dashboardFieldResult.displaySignal;
        if numel(displayed) ~= numel(dashboardFieldResult.tSeconds)
            displayed = real(dashboardFieldResult.rawSignal);
        end
        regular = true(size(dashboardFieldResult.tSeconds));
        regular(dashboardFieldResult.outlierIndices) = false;
        plot(ax, dashboardFieldResult.tSeconds(regular), displayed(regular), 'o', ...
            'DisplayName', 'used');
        if any(~regular)
            plot(ax, dashboardFieldResult.tSeconds(~regular), displayed(~regular), 'x', ...
                'MarkerSize', 8, 'LineWidth', 1.3, 'DisplayName', 'downweighted');
        end
        if isfinite(dashboardFieldResult.R1) && ...
                isfinite(dashboardFieldResult.B) && ...
                isfinite(dashboardFieldResult.C)
            denseTime = linspace(min(dashboardFieldResult.tSeconds), ...
                max(dashboardFieldResult.tSeconds), 300).';
            denseFit = dashboardFieldResult.B - dashboardFieldResult.C .* ...
                exp(-dashboardFieldResult.R1 .* denseTime);
            plot(ax, denseTime, denseFit, '-', 'LineWidth', 1.2, 'DisplayName', 'fit');
        end
        yline(ax, 0, '--');
        grid(ax, 'on');
        xlabel(ax, 'Evolution time (s)');
        ylabel(ax, 'Signed signal (AU)');
        title(ax, sprintf('%d: %.4g mT | T1 %.4g ms | %s | %s', ...
            dashboardFieldResult.fieldIndex, dashboardFieldResult.field_mT, ...
            dashboardFieldResult.T1ms, shortMethod(dashboardFieldResult.method), ...
            shortQc(dashboardFieldResult.qc)), ...
            'Interpreter', 'none');
        hold(ax, 'off');
    end
    title(layout, ['Multi-field recovery fits; crosses are robustly ', ...
        'downweighted points']);
end

function value = shortMethod(method)
    if strcmp(method, 'complex_affine_variable_projection')
        value = 'complex';
    elseif strcmp(method, 'magnitude_all_zero_crossings')
        value = 'magnitude constrained';
    elseif strcmp(method, 'complex_phase_drift_ir')
        value = 'phase-corrected IR';
    elseif strcmp(method, 'magnitude_rician_ir')
        value = 'magnitude IR';
    elseif strcmp(method, 'joint_multifield_inversion_constrained')
        value = 'joint inversion';
    else
        value = method;
    end
end

function value = shortQc(qc)
    if strncmp(qc, 'good', 4)
        value = 'good';
    elseif strncmp(qc, 'failed', 6)
        value = 'failed';
    else
        value = 'review';
    end
end


% ========================================================================
% Small utilities and result templates
% ========================================================================

function r = emptyFieldResult()
    r = struct( ...
        'fieldIndex', NaN, ...
        'field_mT', NaN, ...
        'ok', false, ...
        'method', '', ...
        'tSeconds', [], ...
        'rawSignal', [], ...
        'displaySignal', [], ...
        'displayFit', [], ...
        'predictionComplex', [], ...
        'residual', [], ...
        'weights', [], ...
        'outlierIndices', [], ...
        'B', NaN, ...
        'C', NaN, ...
        'S0', NaN, ...
        'Sinf', NaN, ...
        'R1', NaN, ...
        'T1ms', NaN, ...
        'R1CI95', [NaN NaN], ...
        'T1CI95', [NaN NaN], ...
        'R1SE', NaN, ...
        'T1SE', NaN, ...
        'betaComplex', [NaN NaN], ...
        'phaseRad', NaN, ...
        'phaseTrendRad', [], ...
        'phaseResidualDegrees', NaN, ...
        'orthogonalOffset', NaN, ...
        'orthogonalNRMSE', NaN, ...
        'noiseSigma', NaN, ...
        'dynamicSNR', NaN, ...
        'sse', NaN, ...
        'rsquare', NaN, ...
        'rmse', NaN, ...
        'nrmse', NaN, ...
        'dfe', NaN, ...
        'boundHit', false, ...
        'ciOpen', [false false], ...
        'timeCoverageT1', [NaN NaN], ...
        'zeroCrossingSeconds', NaN, ...
        'signChangeIndex', NaN, ...
        'reliableSignTransitions', NaN, ...
        'signAmbiguous', false, ...
        'candidateRates', [], ...
        'candidateScores', [], ...
        'magnitudePrediction', [], ...
        'magnitudeNoiseModel', '', ...
        'magnitudeNoiseSigma', NaN, ...
        'magnitudeNoiseFloor', NaN, ...
        'magnitudeCoilCount', NaN, ...
        'noiseSource', '', ...
        'selectionReason', '', ...
        'complexCandidateSummary', struct(), ...
        'magnitudeCandidateSummary', struct(), ...
        'independentFitSummary', struct(), ...
        'jointSharedM0', NaN, ...
        'jointFieldRatio', NaN, ...
        'jointAlpha', NaN, ...
        'jointTimeWeights', [], ...
        'jointTotalWeights', [], ...
        'qc', 'not fitted', ...
        'qcReasons', {{}}, ...
        'roiMeta', struct());
end

function c = emptyMagnitudeCandidate()
    c = struct('R1', NaN, 'beta', [NaN; NaN], 'ySigned', [], ...
        'prediction', [], 'weights', [], 'score', Inf, ...
        'boundHit', false, 'signChangeIndex', NaN);
end

function value = safeRsquare(sse, sst)
    if isfinite(sse) && isfinite(sst) && sst > 0
        value = 1 - sse / sst;
    else
        value = NaN;
    end
end

function T1CI = rateCiToT1Ci(rateCI)
    if numel(rateCI) ~= 2 || any(~isfinite(rateCI)) || any(rateCI <= 0)
        T1CI = [NaN NaN];
    else
        T1CI = sort(1000 ./ rateCI);
    end
end

function standardError = intervalToSE(interval, zCritical)
    if numel(interval) == 2 && all(isfinite(interval)) && ...
            isfinite(zCritical) && zCritical > 0
        standardError = (max(interval) - min(interval)) / (2 * zCritical);
    else
        standardError = NaN;
    end
end

function tZero = zeroCrossingTime(B, C, R1)
    tZero = NaN;
    if ~isfinite(B) || ~isfinite(C) || ~isfinite(R1) || C <= 0 || R1 <= 0
        return
    end
    ratio = B / C;
    if ratio > 0 && ratio <= 1
        tZero = -log(ratio) / R1;
    end
end

function tf = memberExists(value, name)
    if isobject(value)
        tf = isprop(value, name);
    else
        tf = isstruct(value) && isfield(value, name);
    end
end

function value = safeSet(value, name, newValue)
    try
        if isobject(value)
            if isprop(value, name)
                value.(name) = newValue;
            end
        elseif isstruct(value)
            value.(name) = newValue;
        end
    catch setException
        warning('T1dispersion:ResultNotStored', ...
            'Could not store result member "%s": %s', name, setException.message);
    end
end

end








    end


end
