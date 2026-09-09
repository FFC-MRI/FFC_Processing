classdef H9_se_nav_v10_cardiac < PulseSequence
    %H9_ir_se Inversion recovery pulse sequence
    
    properties
        typicalSignalLevel  % typical level of signal from the navigator data, to check for signal loss
        initialFactor    % initial offset at the start of an image acquisition, for reversal of corrections
        secondEcho      % set to 1 if using a navigator echo
    end
    
    methods
        function pulse = H9_se_nav_v10_cardiac
            pulse@PulseSequence;
            pulse.pulseSequenceName = 'Inversion Recovery Spin Echo with Triggering';
            pulse.pprFile = 'C:\Scanner\H9_scanner\PulseSequences\H9_se_nav_v10_cardiac.PPR';
            pulse.dBdt = 7;
            pulse.useMode = 3;
            pulse.steppedAcquisition = 0;
            pulse.waveformProfile.pulse = pulse;
            pulse.waveformProfile.generator = InversionRecoveryGenerator;
            pulse.waveformProfile.generator.Tevo = 0;
            pulse.waveformProfile.generator.Bevo = 0.1;
            pulse.waveformProfile.generator.T1est = 0;
            pulse.waveformProfile.generator.indexList = [1,1];
            pulse.waveformProfile.generator.waveform = pulse.waveformProfile;
            pulse.pprParamList = GetPprData(pulse.pprFile);
            pulse.typicalSignalLevel = 0;
            pulse.initialFactor = 0;
            pulse.preAquisitionDelay = 10e-3;
            pulse.secondEcho = 1;
        end
        
        function out = PreparePulseSequence(pulse)
            if isempty(pulse.scan.noiseSample)
                disp('Error: This sequence requires a sample of noise measured beforehand. Run AquireNoise once to do this.')
                out = 0;
                return
            end
            switch pulse.scan.acquisitionMode
                case 'Acquiring'  % acquiring
                    pulse.data = [];    % prepare a slot for the data acquisition
                case 'Setup'
            end
            if isempty(pulse.scan.displayWindow.Children)
                pulse.scan.displayWindow.Children = axes;  % prepare the axes in the main figure for data display
            end
            out = 1;
        end
        
        % this function generates the waveform profile. It should be
        % overriden by the classes derived from this one if needs be.
        function out = MakeWaveformProfile(pulse)
            pulse.waveformProfile.waveList = {};
            UpdatePprParameterList(pulse);
            % readout parameters
            acquisitionField = GetFrequencyInHz(pulse.pprParamList)/gammaH;
            echoTime = GetPprParameter(pulse.pprParamList,'te',0)*1e-3;
            npts = GetPprParameter(pulse.pprParamList,'no_samples',64);
            period = GetPeriodInSecond(pulse.pprParamList);
            tramp =  GetPprParameter(pulse.pprParamList,'tramp',1000)*1e-6;
            t180deg = (4000*1e-6);
            acquisitionTime = (npts*period + echoTime + npts*period/2 + tramp*2+2000e-6+t180deg)*1.1;
            
            pulse.pprParamList = SetPprData(pulse.pprParamList,'no_times',size(pulse.waveformProfile.generator.Tevo,2));
            
            nView = GetPprParameter(pulse.pprParamList,'no_views',1);
            %             iterationNumber =  GetPprParameter(pulse.pprParamList,'effective_no_views',nView); % make sure we undersample if needs be
            iterationNumber =  nView;
            
            pulse.waveformProfile.generator.prePulseDelay = pulse.preAquisitionDelay;
            pulse.waveformProfile.generator.Bread = acquisitionField;
            pulse.waveformProfile.generator.Tread = acquisitionTime;
            pulse.waveformProfile.generator.iterationNumber = iterationNumber;
            pulse.waveformProfile.generator.Tinv = t180deg;
            pulse.waveformProfile.generator.averageNumber = GetPprParameter(pulse.pprParamList,'no_averages',1);
            
            [okflag,pulse.waveformProfile.generator] = Update(pulse.waveformProfile.generator);
            
            
            
            if ~okflag
                error('Invalid waveform. Cannot proceed with the sequence.')
            end
            out = 1;
        end
        
        function out = AcquisitionFunction(pulse,data)
            % This sequence acquires two echoes per acquisition: the first
            % one is an FID and can be used to calibrate the field, the
            % second one is the image data. Each echo triggers this
            % function once.
            out = AcquisitionFunction@PulseSequence(pulse,data);
            no_fields = size(pulse.waveformProfile.generator.Bevo,1);
            no_times = size(pulse.waveformProfile.generator.Tevo,2);
            echonumber = data.Header{5}+1;
            if (echonumber==1)
                linenumber = data.Header{2}+1;
                experimentnumber = double(data.Header{6}'+1);
                current_acquisition = double(data.Header{9}'+1)./2;
                current_time = mod(floor((current_acquisition-1)/GetPprParameter(pulse.pprParamList,'no_views',64)),no_times)+1;
                current_field = mod(floor((current_acquisition-1)/(GetPprParameter(pulse.pprParamList,'no_views',64)*no_times)),no_times)+1;
                total_acquisitions = 2*no_times*no_fields*GetPprParameter(pulse.pprParamList,'no_views',64);
            end
            
            if ((double(data.Header{9}'+1)/2) ==1)
                pulse.userData.start_time=[];
                pulse.userData.start_time=clock;
                tic
                pulse.userData.interval =[];
                pulse.userData.dummymat = zeros(GetPprParameter(pulse.pprParamList,'no_samples',64),GetPprParameter(pulse.pprParamList,'no_views',64),no_times,no_fields);
            end
            
            if ~isfield(pulse.userData,'dummymat')
                pulse.userData.dummymat = zeros(GetPprParameter(pulse.pprParamList,'no_samples',64),GetPprParameter(pulse.pprParamList,'no_views',64),no_times,no_fields);
            end
            
            if (echonumber==1)
                pulse.userData.interval = ([pulse.userData.interval toc]);
                total_time_remaining =  median(pulse.userData.interval).*(total_acquisitions/2-current_acquisition);
                estimated_end = datetime('now') + seconds(total_time_remaining);
                pulse.userData.dummymat(:,linenumber,current_time,current_field) = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
                
                pulse.scan.displayWindow.Children;
                
                subplot(1,2,1)
                imagesc(log(abs(pulse.userData.dummymat(:,:,current_time,current_field)))); axis square off; colormap('gray');
                subplot(1,2,2)
                
                imagesc(rot90(abs(fft2c( pulse.userData.dummymat(:,:,current_time,current_field))),1));  axis square off; colormap('gray');
                delete(findall(gcf,'Tag','stream'));
                dim1 = GetPprParameter(pulse.pprParamList,'no_samples',64);
                dim2 = GetPprParameter(pulse.pprParamList,'no_views',64);
                
                str = {['Evolution Field = ' num2str(pulse.waveformProfile.generator.Bevo(current_field)), ' T '],['Evolution Field = ' num2str(current_field), ' of ', num2str(no_fields)],['Evolution Time = ' num2str(current_time), ' of ', num2str(no_times)],...
                    ['Time remaining = ' num2str(round(total_time_remaining)), ' seconds'],['End time: ' datestr(estimated_end,'HH:MM')] };
                text(-1.2*dim1,1.4*dim2,1,str,'Tag','stream');
                drawnow;
                tic

            end
            if pulse.secondEcho % frequency offset correction if the navigator echo is setup
                if (data.Header{8}==0)&&(data.Header{5}==1)  % only correct over the data from receiver 1 and first echo
                    if GetPprParameter(pulse.pprParamList,'offsetcorr',0)  % check  if the temperature correction slider is set
                        % correct the frequency offset
                        rawData = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
                        % generate the frequency axis
                        T = GetPeriodInSecond(pulse.pprParamList);
                        %                 this should include the FID acquired from all the
                        %                 channels, and should be done once all the channels have
                        %                 been transmitted
                        [offsetInHz,~,successFlag] = ...
                            FindOffsetFromFid(rawData,T,pulse.scan.noiseSample);
                        if successFlag
                            disp(['Frequency offset: ' num2str(offsetInHz) ' Hz.'])
                            if abs(offsetInHz) > 120
                                switch  SlowTemperatureCorrection(pulse.scan,offsetInHz,0.2)
                                    case 0
                                        disp('Error during the frequency offset correction.');
                                        AbortEvo(pulse.scan);
                                    case 1
                                    case 2
                                end
                            end
                        else
                            disp(['Signal too low for correction. Frequency offset estimated: ' num2str(offsetInHz) ' Hz.'])
                        end
                    end
                elseif (data.Header{8}==0)&&(data.Header{5}==size(pulse.scan.noiseSample,9)) % plot the data
                    rawData = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
                    T = GetPeriodInSecond(pulse.pprParamList);
                    %                     subplot(1,2,1);
                    plot(pulse.scan.displayWindow.Children,(1:length(rawData))*T,real(rawData));
                    title(pulse.scan.displayWindow.Children,['Line ' num2str(data.Header{2}+1)])
                    %                     subplot(1,2,2);
                    %                     imagesc(abs(fftshift(ifft2c(pulse.data(:,:,1,1,2,5,1,1,2))))); axis off square; colormap('gray');
                    hold off
                end
            end
            out = 1&out;
        end
        
        function out = FinishAcquisition(pulse)
            if ~isempty(pulse.scan)
                mode = pulse.scan.acquisitionMode;
            else
                mode = 'Acquiring';
            end
            if isequal(mode,'Acquiring')
                data = single(pulse.data); % saving space by avoiding unecessary precision
                %                 try
                %                     pulse.dataProcessed = iterative_images_correction_v8(data,0.15,50,1e-6,[],0);
                %                     figure('WindowStyle','Docked','Visible','on')
                %                     imshow(squeeze(abs(pulse.dataProcessed(:,:,1,1,1,1,1,1,1)))',[0 max(abs(pulse.dataProcessed(:)))])
                %                     disp(['Standard deviation: ' num2str(std(pulse.dataProcessed(:)))])
                %                 end
            end
            % fill the k-space with 0 if compressed acquisitions were used
            linesAcquired = pulse.waveformProfile.generator.kspacelineaquired;
            data = zeros(GetPprParameter(pulse.pprParamList,'no_samples',64),...
                GetPprParameter(pulse.pprParamList,'effective_no_views',64),...
                size(pulse.data,3),...
                size(pulse.data,4),...
                size(pulse.data,5),...
                size(pulse.data,6),...
                size(pulse.data,7),...
                size(pulse.data,8));
            try
                data(:,linesAcquired,:) = pulse.data(:,:,:);
                pulse.data = data;
            catch ME
                disp('Data mismatch, the sequence has probably been aborted.')
            end
            out = 1;
        end
        
        function [out,pulse] = ProcessData(pulse)
            %             pulse.dataProcessed = iterative_images_correction_v8(pulse.data(:,:,:,:,1,:,:),0.15,50,1e-6,[],0);
            out = 1;
        end
        
        function out = ShowData(pulse)
            
            out = 1;
        end
    end
    
    methods (Static)
        
        function file = Reorder(file)
            
            temp1 = load(file.filelocation);
            matfile = temp1.saveList{file.fileindex};
            
            if isempty(matfile.dataProcessed)||(file.reprocess==1)
                
                
                A = double(file.rawdata);
                for p =2:size(A,3)
                    for n=1:size(A,2)
                        thresh = std(A(:,:,:));
                        thresh = mean(squeeze(thresh(:,n,1:2:end)));
                        if std(A(:,n,p))>thresh*5
                            if p ==1
                                A(:,n,p) =  A(:,n,p+1);
                            else
                                A(:,n,p) =  A(:,n,p-1);
                            end
                        end
                    end
                end
                %             file.n_fieldpoints =7; %%%%%%%%%%%%%%%
                A = reshape(A,[file.samples,file.views,file.echoes,file.slices,file.n_timepoints,file.n_fieldpoints]);
                
                if file.echoes>1
                    %                 AA = squeeze(A(:,:,2,:,:,:));
                    
                    AA =A;
                    AA(:,:,2,:,:,:) = [];
                else
                    AA = squeeze(A);
                end
                AA = reshape(AA,[file.samples,file.views,file.slices,file.n_timepoints,file.n_fieldpoints]);
                if file.samples ~= file.views && file.partialkspace ==1
                    AA = padarray(double(AA),[0,double(file.samples-file.views)],0,'post');
                    %             AA = centre_kspace(AA);
                    file.views = file.samples;
                    AA = reshape(AA,file.samples,file.views,[]);
                    
                    %PF recon if we need to
                    try
                        for n=1:size(AA,3)
                            [~, AA(:,:,n)] = pocs(AA(:,:,n),10,0);
                        end
                    catch
                    end
                end
                
                %             file.views = file.samples;
                AA = reshape(AA,[file.samples,file.views,file.slices,file.n_timepoints,file.n_fieldpoints]);
                
                %             h= figure;
                %             imagesc(abs(ifft2c(AA(:,:,1,1,1))));
                %             mask = roipoly;
                %             close(h);
                if file.backgroundselect ==1
                    hh = figure;
                    imagesc(abs(ifft2c(AA(:,:,1,1,1)))); axis off square; colormap('gray');
                    background = roipoly;
                    close(hh);
                else
                    background = [];
                end
                corrected_stack = (iterative_images_correction_v7(AA,0.20,20,0.00015,background));
                %             corrected_stack = register_images(corrected_stack);
                if file.samples ~= file.views
                    corrected_stack = padarray(double(corrected_stack),[0,double(file.samples-file.views)/2],0,'both');
                    file.views = file.views;
                end
                file.complexkspace = fft2c(corrected_stack);
                file.originalcomplexkspace =file.complexkspace;
                matfile.dataProcessed = file.complexkspace;
                temp1.saveList{file.fileindex} = matfile;
                saveList = temp1.saveList;
                [directory,~] = bst_fileparts(file.filelocation);
                processed_dir = fullfile(directory,'ProcessedData.mat');
                save(processed_dir,'saveList');
            else
                file.complexkspace = matfile.dataProcessed;
                file.originalcomplexkspace =file.complexkspace;
            end
        end
    end
end

