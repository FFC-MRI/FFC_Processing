classdef H9_ir_se_nav_v4 < PulseSequence
    %H9_ir_se Inversion recovery pulse sequence
    
    properties
        typicalSignalLevel  % typical level of signal from the navigator data, to check for signal loss
        initialFactor    % initial offset at the start of an image acquisition, for reversal of corrections
        secondEcho      % set to 1 if using a navigator echo
    end
    
    methods
        function pulse = H9_ir_se_nav_v4
            pulse@PulseSequence;
            pulse.pulseSequenceName = 'Inversion Recovery Spin Echo';
            pulse.pprFile = 'C:\Scanner\H9_scanner\PulseSequences\H9_se_nav_v4.PPR';
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
            t180deg = (1500*1e-6)+GetPprParameter(pulse.pprParamList,'p90',120)*3e-6; %180 is 3x90 duration
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
            pulse.dataProcessed = iterative_images_correction_v8(pulse.data(:,:,:,:,1,:,:),0.15,50,1e-6,[],0);
            out = 1;
        end
        
        function out = ShowData(pulse)
            if isobject(pulse.scan)
                figure(pulse.scan.displayWindow)
                set(pulse.scan.displayWindow,'Visible','on')
            else
                figure
            end
            nimage = size(pulse.dataProcessed,7);
            nline = sqrt(nimage);
            if ceil(nline)*floor(nline) >= nimage
                n2 = floor(nline);
            end
            nline = ceil(nline);
            for indimage = 1:nimage
                subplot(n2,nline,indimage);
                imshow(squeeze(abs(pulse.dataProcessed(:,:,1,1,1,1,indimage,1,1)))',[0 max(abs(pulse.dataProcessed(:)))])
            end
            truesize(gcf,[600 600]);
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
                AA(:,:,1,:,:,:) = [];
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
                hh = figure;
                imagesc(abs(ifft2c(AA(:,:,1,1,1)))); axis off square; colormap('gray');
                background = roipoly;
                close(hh);
            else
                background = [];
            end
            corrected_stack = (iterative_images_correction_v8(AA,0.20,50,0.00015,background,0));
            stackCorrectedFft = fftshift(fftshift(fft(fft(corrected_stack,[],1),[],2),1),2);
            corrected_stack = ifft(ifft(stackCorrectedFft,[],1),[],2);
%             corrected_stack = register_images(corrected_stack);
            
            file.complexkspace = fft2c(corrected_stack);
            file.originalcomplexkspace =file.complexkspace;
            matfile.dataProcessed = file.complexkspace;
            temp1.saveList{file.fileindex} = matfile;
            saveList = temp1.saveList;
            save('PulseSequenceList.mat','saveList');
            else
                file.complexkspace = matfile.dataProcessed;
                file.originalcomplexkspace =file.complexkspace;
            end
        end
    end
end

