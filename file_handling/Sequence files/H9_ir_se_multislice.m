classdef H9_ir_multislice < PulseSequence
    %H9_ir_se Inversion recovery pulse sequence
    
    properties
        %         targetFrequency     % frequency of the first line of an image, used for temperature correction
        typicalSignalLevel  % typical level of signal from the navigator data, to check for signal loss
        initialFactor    % initial offset at the start of an image acquisition, for reversal of corrections
        secondEcho      % set to 1 if using a navigator echo
    end
    
    methods
        function pulse = H9_ir_se_multislice
            pulse@PulseSequence;
            pulse.pulseSequenceName = 'Inversion Recovery Multislice spin echo';
            pulse.pprFile = 'C:\Scanner\H9_scanner\PulseSequences\H9_ir_se_multislice.PPR';
            pulse.dBdt = 7;
            pulse.useMode = 3;
            pulse.steppedAcquisition = 0;
            pulse.waveformProfile.pulse = pulse;
            pulse.waveformProfile.generator.Tevo = 0;
            pulse.waveformProfile.generator.Bevo = 0.1;
            pulse.waveformProfile.generator.T1est = 0.1;
            pulse.waveformProfile.generator.indexList = [1,1];
            pulse.waveformProfile.generator.inversionFlag = false;
            pulse.pprParamList = GetPprData(pulse.pprFile);
            %             pulse.targetFrequency = 0;
            pulse.typicalSignalLevel = 0;
            pulse.initialFactor = 0;
            pulse.preAquisitionDelay = 10e-3;
            pulse.secondEcho = 0;
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
            out = UpdatePprParameterList(pulse);
            % Calculation of delays
            % inversion parameters
            if GetPprParameter(pulse.pprParamList,'hard_inv',0)
                inversionPulseDuration = GetPprParameter(pulse.pprParamList,'p180',0)*1e-6;
            else
                switch GetPprParameter(pulse.pprParamList,'rfnum_inv',5)
                    case 1
                        inversionPulseDuration = 1332e-6;
                    case 2
                        inversionPulseDuration = 2664e-6;
                    case 3
                        inversionPulseDuration = 5328e-6;
                    case 4
                        inversionPulseDuration = 2000e-6;
                    case 5
                        inversionPulseDuration = 4000e-6;
                    case 6
                        inversionPulseDuration = 8000e-6;
                    case 7
                        inversionPulseDuration = 2000e-6;
                    case 8
                        inversionPulseDuration = 4000e-6;
                    case 9
                        inversionPulseDuration = 2000e-6;
                                            
                end
            end
            gradRampTime =  GetPprParameter(pulse.pprParamList,'tramp',1000)*1e-6*1.1;  % adding 10 % time for good measure
            inversionPulseDuration = inversionPulseDuration + gradRampTime*2;
            % pulse durations
            switch GetPprParameter(pulse.pprParamList,'rfnum',5)
                case 1
                    pulseDuration = 1332e-6;
                case 2
                    pulseDuration = 2664e-6;
                case 3
                    pulseDuration = 5328e-6;
                case 4
                    pulseDuration = 2000e-6;
                case 5
                    pulseDuration = 4000e-6;
                case 6
                    pulseDuration = 8000e-6;
                case 7
                    pulseDuration = 2000e-6;
                case 8
                    pulseDuration = 4000e-6;
                case 9
                    pulseDuration = 2000e-6;

            end
            % The prepulse delay is already taken into accound within the
            % WaveformProfile object.
            % readout parameters
            acquisitionField = GetFrequencyInHz(pulse.pprParamList)/gammaH;
            npts = GetPprParameter(pulse.pprParamList,'no_samples',64);
            period = GetPeriodInSecond(pulse.pprParamList);
            trdd = GetPprParameter(pulse.pprParamList,'trdd',0)*1e-6;
            postpulseDelay = GetPprParameter(pulse.pprParamList,'taftfc2',0)*1e-3;
            echoTime = GetPprParameter(pulse.pprParamList,'te',0)*1e-3;
            echoNumber = GetPprParameter(pulse.pprParamList,'no_echoes',1);
            sliceNumber = GetPprParameter(pulse.pprParamList,'no_slices',1);
            acquisitionTime = (pulse.preAquisitionDelay + postpulseDelay + trdd +  (echoTime + pulseDuration/2 + npts*period + gradRampTime*2)*sliceNumber)*1.1; % add 10% time for good measure
            %             % Detect if an evolution period is required:
            %             evolutionSwitch = ~isempty(pulse.waveformParam.evolutionField) ...
            %                               && pulse.waveformParam.evolutionTime(1)>0;
            
            nView = GetPprParameter(pulse.pprParamList,'no_views',1);
            iterationNumber = nView;
            % crop the settling time before the acquisition by merging it
            % with the EVO wainting time
            pulse.waveformProfile.generator.prePulseDelay = pulse.preAquisitionDelay;
            
            pulse.waveformProfile.generator.Bread = acquisitionField;
            pulse.waveformProfile.generator.Tread = acquisitionTime;
            pulse.waveformProfile.generator.inversionTime = inversionPulseDuration;
            pulse.waveformProfile.experimentNumber = iterationNumber;
            pulse.waveformProfile.generator.iterationNumber = iterationNumber;
            % Detect if an evolution period is required:
            if GetPprParameter(pulse.pprParamList,'inv_obs_mod_level',0); 
                pulse.waveformProfile.generator.inversionFlag = true;
            else
                pulse.waveformProfile.generator.inversionFlag = false;
            end
                              
            % finally, sort out the number of averages
            pulse.waveformProfile.generator.averageNumber = GetPprParameter(pulse.pprParamList,'no_averages',1);
            
            [okflag,pulse.waveformProfile.generator] = Update(pulse.waveformProfile.generator);
            if ~okflag
                error('Invalid waveform. Cannot proceed with the sequence.')
            end            
        end
        
        function out = AcquisitionFunction(pulse,data)
            % This sequence acquires two echoes per acquisition: the first
            % one is an FID and can be used to calibrate the field, the
            % second one is the image data. Each echo triggers this
            % function once.
            out = AcquisitionFunction@PulseSequence(pulse,data);
            %             if pulse.secondEcho % frequency offset correction if the navigator echo is setup
            %                 if (data.Header{8}==0)&&(data.Header{5}==0)  % only correct over the data from receiver 1 and first echo
            %                     if GetPprParameter(pulse.pprParamList,'offsetcorr',0)  % check  if the temperature correction slider is set
            %                         % correct the frequency offset
            %                         rawData = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
            %                         % generate the frequency axis
            %                         T = GetPeriodInSecond(pulse.pprParamList);
            %                 %                 this should include the FID acquired from all the
            %                 %                 channels, and should be done once all the channels have
            %                 %                 been transmitted
            %                         [offsetInHz,~,successFlag] = ...
            %                         FindOffsetFromFid(rawData,T,pulse.scan.noiseSample);
            %                         if successFlag
            %                             disp(['Frequency offset: ' num2str(offsetInHz) ' Hz.'])
            %                             if abs(offsetInHz) > 120
            %                                 switch  SlowTemperatureCorrection(pulse.scan,offsetInHz,0.2)
            %                                     case 0
            %                                         disp('Error during the frequency offset correction.');
            %                                         AbortEvo(pulse.scan);
            %                                     case 1
            %                                     case 2
            %                                 end
            %                             end
            %                         else
            %                             disp(['Signal too low for correction. Frequency offset estimated: ' num2str(offsetInHz) ' Hz.'])
            %                         end
            %                     end
            %                 elseif (data.Header{8}==0)&&(data.Header{5}==size(pulse.scan.noiseSample,9)) % plot the data
            %                     rawData = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
            %                     T = GetPeriodInSecond(pulse.pprParamList);
            % %                     subplot(1,2,1);
            %                     plot(pulse.scan.displayWindow.Children,(1:length(rawData))*T,real(rawData));
            %                     title(pulse.scan.displayWindow.Children,['Line ' num2str(data.Header{2}+1)])
            % %                     subplot(1,2,2);
            % %                     imagesc(abs(fftshift(ifft2c(pulse.data(:,:,1,1,2,5,1,1,2))))); axis off square; colormap('gray');
            %                     hold off
            %                 end
            %             end
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
                try
                    pulse.dataProcessed = iterative_images_correction_v8(data,0.15,50,1e-6,[],0);
                    figure('WindowStyle','Docked','Visible','on')
                    imshow(squeeze(abs(pulse.dataProcessed(:,:,1,1,1,1,1,1,1)))',[0 max(abs(pulse.dataProcessed(:)))])
                    disp(['Standard deviation: ' num2str(std(pulse.dataProcessed(:)))])
                end
            end
            out = 1;
        end
        
    end
    
   methods (Static) 
    
          function file = Reorder(file)
            
            temp1 = load(file.filelocation);
            matfile = temp1.saveList{file.fileindex};
            
            if isempty(matfile.dataProcessed)||(file.reprocess==1)
            
            
            A = double(file.rawdata);
           for p =1:size(A,3)
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
            if file.samples ~= file.views
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
            
            file.views = file.samples;
            AA = reshape(AA,[file.samples,file.views,file.slices,file.n_timepoints,file.n_fieldpoints]);
            
            %             h= figure;
            %             imagesc(abs(ifft2c(AA(:,:,1,1,1))));
            %             mask = roipoly;
            %             close(h);
            if file.backgroundselect ==1
                for n=1:file.slices
                hh = figure;
                imagesc(abs(ifft2c(AA(:,:,n,1,1)))); axis off square; colormap('gray');
                background(:,:,n) = roipoly;
                close(hh);
                end
            else
                background = [];
            end
            
            for l=1:file.slices    
            corrected_stack(:,:,l) = iterative_images_correction_v7(AA(:,:,l),0.25,20,1e-6,[]);
            end
            
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
                file.complexkspace = squeeze(matfile.dataProcessed);
                file.originalcomplexkspace =squeeze(file.complexkspace);
            end
        end
    end
end
            

