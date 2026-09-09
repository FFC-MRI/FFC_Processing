classdef h9_flash_v2 < PulseSequence
    %H9_ir_se Inversion recovery pulse sequence
    
    properties
        
    end
    
    methods
        function pulse = h9_flash_v2
            pulse@PulseSequence;
            pulse.pulseSequenceName = 'Turbo Flash';
            pulse.pprFile = 'C:\Scanner\H9_scanner\PulseSequences\h9_flash_v2.PPR';
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
            pulse.preAquisitionDelay = 10e-3;
        end
        
        % this function generates the waveform profile. It should be
        % overriden by the classes derived from this one if needs be.
        function out = MakeWaveformProfile(pulse)
            out = UpdatePprParameterList(pulse);
            % reinitialise the waveform object
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
                otherwise
                         pulseDuration = 2000e-6;

            end
            % readout parameters
            polarisationTime = GetPprParameter(pulse.pprParamList,'t_pol',100)*1e-3;
            polarisationField = GetPprParameter(pulse.pprParamList,'b_pol',200)*1e-3;
            acquisitionField = GetFrequencyInHz(pulse.pprParamList)/gammaH;
            npts = GetPprParameter(pulse.pprParamList,'no_samples',64);
            tr =  GetPprParameter(pulse.pprParamList,'tr',1)*1e-3;
            period = GetPeriodInSecond(pulse.pprParamList);
            trdd = GetPprParameter(pulse.pprParamList,'trdd',70)*1e-6;
            postpulseDelay = GetPprParameter(pulse.pprParamList,'taftfc2',0)*1e-3 + 0.05;
            echoTime = GetPprParameter(pulse.pprParamList,'te',0)*1e-3;
            echoNumber = GetPprParameter(pulse.pprParamList,'no_echoes',1);
            lineNumber = GetPprParameter(pulse.pprParamList,'no_views',1);
            tramp = GetPprParameter(pulse.pprParamList,'tramp',1)*1e-6;
            tref = GetPprParameter(pulse.pprParamList,'tref',1)*1e-6;
            tdp = GetPprParameter(pulse.pprParamList,'tdp',1)*1e-6;
            t180deg = GetPprParameter(pulse.pprParamList,'p90',1)*2e-6;
            no_experiments = GetPprParameter(pulse.pprParamList,'no_experiments',1);
            spoilerTime = 800e-6;
            acquisitionTime = (npts*period + trdd + echoTime*echoNumber + tramp*8 + 89e-6 + tref*1.3 + pulseDuration + tdp + spoilerTime)*1.2*lineNumber + postpulseDelay; % add 10% time for good measure
            % recycle delay parameters
            recycleDelay = acquisitionTime*1.5;
            % make the waveforms
            nView2 = GetPprParameter(pulse.pprParamList,'no_views_2',1);
            sliceNumber = GetPprParameter(pulse.pprParamList,'no_slices',1);
            if GetPprParameter(pulse.pprParamList,'inversion_onoff',0)
            pulse.waveformProfile.generator.Tinv =t180deg;
            else
            pulse.waveformProfile.generator.Tinv =0;
            end
            iterationNumber = nView2*sliceNumber;     
            
            pulse.waveformProfile.generator.prePulseDelay = pulse.preAquisitionDelay;
            pulse.waveformProfile.generator.Bread = acquisitionField;
            pulse.waveformProfile.generator.Tread = acquisitionTime;
            pulse.waveformProfile.generator.iterationNumber = iterationNumber;
            pulse.waveformProfile.generator.averageNumber = GetPprParameter(pulse.pprParamList,'no_averages',1);
            
            [okflag,pulse.waveformProfile.generator] = Update(pulse.waveformProfile.generator);
            
 
            if ~okflag
                error('Invalid waveform. Cannot proceed with the sequence.')
            end
            
            out = 1;
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
                    pulse.data = [];  % acquire data anyway
            end
            figure(pulse.scan.displayWindow)
            subplot(3,3,1)
            out = 1;
        end
        
        function out = AcquisitionFunction(pulse,data)
            % This sequence acquires two echoes per acquisition: the first
            % one is an FID and can be used to calibrate the field, the
            % second one is the image data. Each echo triggers this
            % function once.  
            out = AcquisitionFunction@PulseSequence(pulse,data);
            lineNumber = double(data.Header{2}+1);
            sliceNumber = double(data.Header{4}+1);
            if lineNumber == GetPprParameter(pulse.pprParamList,'no_views',1)
                switch pulse.scan.acquisitionMode
                    case 'Setup'
                        table = [0,-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8,8,-9,9,-10,10,-11,11,-12,12,-13,13,-14,14,-15,15,-16,16,-17,17,-18,18,-19,19,-20,20,-21,21,-22,22,-23,23,-24,24,-25,25,-26,26,-27,27,-28,28,-29,29,-30,30,-31,31,-32,32,-33,33,-34,34,-35,35,-36,36,-37,37,-38,38,-39,39,-40,40,-41,41,-42,42,-43,43,-44,44,-45,45,-46,46,-47,47,-48,48,-49,49,-50,50,-51,51,-52,52,-53,53,-54,54,-55,55,-56,56,-57,57,-58,58,-59,59,-60,60,-61,61,-62,62,-63,63,-64];
                       tempdata = pulse.data(:,:,1,sliceNumber,1,1,1,1,1);
                       
                       tempmatrix = zeros(size(tempdata));
                       for n=1:size(tempdata,2)  
                       tempmatrix(:,1+size(tempdata,2)/2+table(n),:,:,:,:,:,:,:) =  tempdata(:,n,:,:,:,:,:,:,:);
                       end
                        figure(pulse.scan.displayWindow)
                        subplot(1,GetPprParameter(pulse.pprParamList,'no_slices',1),sliceNumber)
                        im = rot90(ifft2c(tempmatrix),-1);
                        imagesc(squeeze(abs(im))); colormap('gray'); axis off square
                end
            end
        end
        
        function out = FinishAcquisition(pulse)
            if isequal(pulse.scan.acquisitionMode,'Acquiring')
%                 data = squeeze(pulse.data);
%                 set(scan.displayWindow,'Visible','on')
%                 figure(scan.displayWindow)
%                 for i = 1:size(data,3)
%                     im = fftshift(ifft2(squeeze(data(:,:,i))));
%                     subplot(1,3,i)
%                     imshow(abs(im),[0 max(abs(im(:)))])
%                 end
%                 drawnow
%                 figure
%                 for i = 1:size(data,1)  
%                     [I1,ph,A1,bkgd,jj] = iterative_images_correction_v7(fftshift(squeeze(data(i,:,:))),0.15,20,0.15,[]);
%                     subplot(1,3,i)
%                     imshow(abs(I1),[0 max(abs(I1(:)))])
%                 end
%                 disp(['Standard deviation: ' num2str(std(im(:)))])
            end
            out = 1;
        end
        
        function [out,pulse] = ProcessData(pulse)
%             pulse.dataProcessed = iterative_images_correction_v8(pulse.data,0.15,50,1e-6,[],0);
            out = 1;
        end
        
        function out = ShowData(pulse)
%             if isobject(pulse.scan)
%                 figure(pulse.scan.displayWindow)
%                 set(pulse.scan.displayWindow,'Visible','on')
%             else
%                 figure
%             end
%             nimage = size(pulse.dataProcessed,4);
%             for indimage = 1:nimage
%                 subplot(1,nimage,indimage);
%                 imshow(squeeze(abs(pulse.dataProcessed(:,:,1,indimage,1,1,1,1,1)))',[0 max(abs(pulse.dataProcessed(:)))])
%             end
%             truesize(gcf,[600 600])
            out = 1;
        end
        
    end
        methods (Static) 
            
        function file = Reorder(file)
            temp1 = load(file.filelocation);
            matfile = temp1.saveList{file.fileindex};
            
            if isempty(matfile.dataProcessed)||(file.reprocess==1)
            
         A = double(file.rawdata); 

%                 for p =1:size(A,3)
%                     for n=1:size(A,2)
%                         thresh = std(A(:,:,:));
%                         thresh = mean(squeeze(thresh(:,n,1:2:end)));
%                         if std(A(:,n,p))>thresh*5
%                             if p ==1
%                                 A(:,n,p) =  A(:,n,p+1);
%                             else
%                                 A(:,n,p) =  A(:,n,p-1);
%                             end
%                         end
%                     end
%                 end
            
            A = reshape(A,[file.samples,file.views,file.echoes,file.slices,file.n_timepoints,file.n_fieldpoints]);  
            
            AA =A;
            
            table = [0,-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8,8,-9,9,-10,10,-11,11,-12,12,-13,13,-14,14,-15,15,-16,16,-17,17,-18,18,-19,19,-20,20,-21,21,-22,22,-23,23,-24,24,-25,25,-26,26,-27,27,-28,28,-29,29,-30,30,-31,31,-32,32,-33,33,-34,34,-35,35,-36,36,-37,37,-38,38,-39,39,-40,40,-41,41,-42,42,-43,43,-44,44,-45,45,-46,46,-47,47,-48,48,-49,49,-50,50,-51,51,-52,52,-53,53,-54,54,-55,55,-56,56,-57,57,-58,58,-59,59,-60,60,-61,61,-62,62,-63,63,-64];
                       tempdata = AA;
                       
                       tempmatrix = zeros(size(tempdata));
                       for n=1:size(tempdata,2)  
                       tempmatrix(:,1+size(tempdata,2)/2+table(n),:,:,:,:,:,:,:) =  tempdata(:,n,:,:,:,:,:,:,:);
                       end
                       
            AA = tempmatrix;
            AA = reshape(AA,[file.samples,file.views,file.slices,file.n_timepoints,file.n_fieldpoints]);
            AA = centre_kspace(AA);
            
            if file.backgroundselect ==1
                hh = figure;
                imagesc(abs(ifft2c(AA(:,:,1,1,1)))); axis off square; colormap('gray');
                background = roipoly;
                close(hh);
            else
                background = [];
            end
%             AA = permute(AA,[2,1,3,4,5,6,7]);
%              corrected_stack = iterative_images_correction_v8(AA,0.20,20,0.0015,background,0);
%              stackCorrectedFft = fftshift(fftshift(fft(fft(corrected_stack,[],1),[],2),1),2);
%              corrected_stack = ifft(ifft(stackCorrectedFft,[],1),[],2);
%           corrected_stack = AA;
% %              corrected_stack = register_images(corrected_stack);
%             file.complexkspace = fft2c*AA;
             file.complexkspace =AA;
           
            file.originalcomplexkspace =file.complexkspace;
            matfile.dataProcessed = file.complexkspace;
            temp1.saveList{file.fileindex} = matfile;
            saveList = temp1.saveList;
            save('ProcessedData.mat','saveList');
            else
                file.complexkspace = matfile.dataProcessed;
                file.originalcomplexkspace =file.complexkspace;
            end
        end
    end
end
