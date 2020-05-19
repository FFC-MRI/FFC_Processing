classdef H9_se_propeller < PulseSequence
    %H9_ir_se Inversion recovery pulse sequence
    
    properties
        typicalSignalLevel  % typical level of signal from the navigator data, to check for signal loss
        initialFactor    % initial offset at the start of an image acquisition, for reversal of corrections
        secondEcho      % set to 1 if using a navigator echo
    end
    
    methods
        function pulse = H9_se_propeller
            pulse@PulseSequence;
            pulse.pulseSequenceName = 'Propellor Spin Echo';
            pulse.pprFile = 'C:\Scanner\H9_scanner\PulseSequences\H9_se_propeller.PPR';
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
            
            nView2 = GetPprParameter(pulse.pprParamList,'no_views_2',1);
            nView = GetPprParameter(pulse.pprParamList,'no_views',1);
            %             iterationNumber =  GetPprParameter(pulse.pprParamList,'effective_no_views',nView); % make sure we undersample if needs be
            iterationNumber =  nView*nView2;
            
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
            
            
            
            
            tic
            nb = file.views;   %#blades
            nl = double(file.views2); %#lines per blade
            s = double(file.samples);
            A = double(file.rawdata);
            AA = permute(A,[1 3 2]);
            winfactor=4;
            win_n=s/winfactor;
            win_m=nl/winfactor;
            wind = (triang(win_n)*triang(win_m)');
            bladeangle = file.bladeangles(2)./10;
            wind_pad = padarray(wind, [floor((s-win_n)/2) floor((nl-win_m)/2)], 'replicate','post');
            wind = padarray(wind_pad, [ceil((s-win_n)/2) ceil((nl-win_m)/2)], 'replicate','pre');
            zeropadfactor =1;
            rotationcorrection =1;
            translationcorrection=1;
            phasecorrection=1;
            method = 'TV';
            % AA = padarray(AA,double([s/4 nl/4]));
            %% Phase Correction
            for n=1:nb
                AAFT =ifft2c(padarray(AA(:,:,n),(zeropadfactor-1)*double([s/2 nl/2])));
                AAFTwind = ifft2c(padarray((AA(:,:,n).*wind),(zeropadfactor-1)*double([s/2 nl/2])));
                %  phase1(:,:,n) = (atan((imag(AAFTwind(:,:,n))./real(AAFTwind(:,:,n)))));
                % phase2(:,:,n)=(atan((imag(AAFT(:,:,n))./real(AAFT(:,:,n)))));
                %phase1(:,:,n)=unwrap(angle(AAFTwind(:,:,n)));
                %   phase2(:,:,n)=angle(AAFT(:,:,n));
                %   phaseinfo=(phase1-phase2);
                if phasecorrection ==1
                    AAFTcorrect=AAFT.*exp(-1i*angle(AAFTwind));
                    AAcorrect(:,:,n)=fft2c(updownsample( AAFTcorrect,nl,s,0,0 ));
                else
                    AAcorrect(:,:,n)=fft2c(AAFT);
                end
                %%
                
                
            end
            s=s*1;
            nl=nl*1;
            
            
%             injectdeltar=0;
%             injectdeltac=10;
%             [nr,nc,~]=size(AAcorrect);
%             [u,v] = meshgrid(1:nc,1:nr);
%             exp_phase_ramp=exp(1i*2*pi*((u.*injectdeltar)/nc+(v.*injectdeltac)/nr));
%             AAcorrect(:,:,1) = ((AAcorrect(:,:,1)).*exp_phase_ramp);
           
            
            
%             %% Bulk Rotation Correction
%             truncatedsignal=AAcorrect(s/2-nl/2+1:s/2+nl/2,:,:);
%             x = 1:nl;             % The range of x values.
%             y = 1:nl;             % The range of y values.
%             [X,Y] = meshgrid (x,y);
%             wmat=sqrt((X-nl/2).^2+(Y-nl/2).^2);           
%          
%             for n=1:nb
%                 k = PropTraj(nl,nl,nb,bladeangle,0,n);              
%                 m =  reshape(truncatedsignal(:,:,n).*1,[nl*nl,1]);
%                 m = double(m.*repmat(exp(1i*pi*(k(:,1)+k(:,2))),[1, 1]));
%                 [a,~,~] = Prepare4Recon(m, k, ones(nl,nl), ones(nl,nl));
%                 bladeim(:,:,n)=abs(a);
%             end
%             referenceim=mean(bladeim,3);
% 
%             for n=1:nb
%                 m1=  reshape(truncatedsignal(:,:,n).*1,[nl*nl,1]);
%                 for rotang=1:180;
%                     ang(rotang)=(rotang-90)./2;     
%                     k = PropTraj(nl,nl,nb,bladeangle,ang(rotang),n);
%                     m = double(m1.*repmat(exp(1i*pi*(k(:,1)+k(:,2))),[1, 1]));
%                     [a,~,~] = Prepare4Recon(m, k, ones(nl,nl), ones(nl,nl));
%                     bladeimtemp=abs(a);
%                     corval(n,rotang) =corr2(bladeimtemp,referenceim);
%                 end
% %                 fit = polyfit(ang,corval(n,:),2);
% %                 finegrid=linspace(min(ang),max(ang),100);
% %                 anglefit = polyval(fit,finegrid);
%                 [~, ind]=max(corval(n,:));
%                 anglecorrection(n)=1*rotationcorrection.*ang(ind);
%             end
            %%
            
            
%             %% BULK TRANSLATION CORRECTION
%             x=1:nl;
%             x = 1:nl;             % The range of x values.
%             y = 1:nl;             % The range of y values.
%             [X,Y] = meshgrid (x,y);
%           
%             for n=1:nb
%                 k = PropTraj(nl,nl,nb,bladeangle,anglecorrection(n),n);
%                 m=  reshape(truncatedsignal(:,:,n),[nl*nl,1]);
%                 m = double(m.*repmat(exp(1i*pi*(k(:,1)+k(:,2))),[1, 1]));
%                 [a(:,:,n),~,~] = Prepare4Recon(m, k, ones(nl,nl), ones(nl,nl));
%                 bladedata(:,:,n)=(fft2(a(:,:,n)));
%             end
%             referencedata=mean(bladedata,3);
%             
%             
%             for n=1:nb
%                 test=fftshift(ifft2((bladedata(:,:,n).*conj(referencedata))));    %phase correlation of average data set with individual blades
%                 [ampl,ind]=max(abs(test(:)));
%                 [x,y]=ind2sub(size(test),ind);      %identify maximum amplitude pixel in correlation image
%                 X= x-1:x+1;
%                 Y = y-1:y+1;
%                 %papers authors use a 3 pt parabolic fit, so we will too
%                 yvals = abs(test(X,y))';       %identify the 3 points along each direction about the maximum
%                 xvals = abs(test(x,Y));
%                 p1 = polyfit(X,xvals,2);        %do parabolic fit in each direction
%                 p2 = polyfit(Y,yvals,2);
%                 
%                 finegridx=linspace(x-1,x+1,201);        %work out the shift in sub pixels. technically this is interpolation so be cautious
%                 finegridy=linspace(y-1,y+1,201);
%                 x1 = polyval(p1,finegridx);
%                 y1 = polyval(p2,finegridy);
%                 
%                 [~, ind]=max(x1);
%                 delx = finegridx(ind);
%                 [~, ind]=max(y1);
%                 dely = finegridy(ind);
%                 
%                 deltar=-1*(nl/2+1-delx);        %this is the actual shift in pixels
%                 deltac=-1*(nl/2+1-dely);
%                 
%                 deltar=deltar*translationcorrection;
%                 deltac=deltac*translationcorrection;
%                 
%                 exp_phase_ramp=exp(1i*2*pi*((u.*deltar)/nc+(v.*deltac)/nr));    %apply the translation in fourier space
%                 AAcorrect(:,:,n) = AAcorrect(:,:,n).*exp_phase_ramp;
%                 
%             end
            %%
            
            %% Correlation Weighting
%             truncatedsignal=AAcorrect(s/2-nl/2+1:s/2+nl/2,:,:);
%             truncatedsignalabs = abs(truncatedsignal);
%             
%             for n=1:nb
%                 reducedk(:,:,n) = imrotate(truncatedsignal(:,:,n),(n-1)*bladeangle+anglecorrection(n),'bicubic','crop');
%                 R2(:,:,n)=reducedk(:,:,n);
%                 R2abs(:,:,n)=abs(R2(:,:,n));
%                 R2im(:,:,n)=(ifft2(R2(:,:,n)));
%             end
%             Daabs=mean(R2abs,3);
%             Da=mean(R2,3);
%             
%             for n=1:12
%                 Xn(:,n)=abs(sum(sum((Da.*conj(R2(:,:,n))))));
%                 
%             end
%             for n=1:12
%                 Pn(:,n)=(0.1+0.9*((Xn(:,n)-min(Xn))/(max(Xn)-min(Xn))))^2;
%             end
%             %%
            
            
            
            
            %%  GRID AND RECONSTRUCTION
            clearvars -except s obj nl nb bladeangle anglecorrection AAcorrect method;
          
            m=reshape((AAcorrect),[s*nl*nb,1]);           
            k=PropTraj(s,nl,nb,bladeangle);           
            m = double(m.*repmat(exp(1i*pi*(k(:,1)+k(:,2))),[1, 1]));
            clearvars -except m k s obj method
            [a,A,P] = Prepare4Recon(m, k, ones(s,s), ones(s,s)); colormap('gray')
            
            switch method
                case 'TV'
                    [AAA,~,~] = ReconTV((a), A, 1.0e-3*max(abs(a(:))), a, 30, 6, P);
                    
                case 'wavelet'
                    mxsize=s*[1 1];
                    W = DWT(3*[1,1],mxsize, true, 'haar');
%                     W = BlockDCT(8); % uncomment this if you prefer the 8x8 BlockDCT as a regularization transform
                    alpha = PowerIteration(A, a);
%                     alpha = PowerIterationWav( @(w) W*(A(W'*(w))), W*a ); % use this one to perform FWISTA
                    % the alpha can depends on the MRI setting, not on the data, it can be
                    % precomputed and saved for future reconstructions involving same k-space
                    % trajectory, coil sensitivities and reconstruction matrix.
                    
                    [AAA,~,~] = ReconWavFISTA(a, A, 1.0e-3*max(abs(a(:))), W, alpha, a, 100, true);
                    
            end
            AAA=(AAA);
            [nr,nc,~]=size(AAA);
            [u,v] = meshgrid(1:nc,1:nr);
            exp_phase_ramp=exp(1i*2*pi*((u.*nr/2)/nc+(v.*nr/2)/nr));
            AAA=(fft2(ifftshift(AAA.*exp_phase_ramp))); %undo whacky origin shift to sync kspace and image space origins
            figure,
            imagesc(abs(ifft2c(AAA)));
            toc
            %%
        
                
                
                
     
                
            
            file.complexkspace = AAA;
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

