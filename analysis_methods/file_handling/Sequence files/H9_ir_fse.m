classdef H9_ir_fse < PulseSequence
    %H9_ir_se Inversion recovery pulse sequence
    
    properties
        
    end
    
    methods
        function pulse = H9_ir_fse
            pulse@PulseSequence;
            pulse.pulseSequenceName = 'Inversion recovery fast spin echo';
            pulse.pprFile = 'C:\Scanner\H9_scanner\PulseSequences\H9_ir_fse.PPR';
            pulse.dBdt = 7;
            pulse.useMode = 3;
            pulse.steppedAcquisition = 0;
            pulse.waveformParam.evolutionTime = 0.1;
            pulse.waveformParam.evolutionField = 0.1;
            pulse.waveformParam.estimatedT1 = 0.1;
            pulse.waveformParam.indexList = [1,1];
            pulse.pprParamList = GetPprData(pulse.pprFile);
        end
        
         function out = MakeWaveformProfile(pulse)
            % readout parameters
            polarisationTime = GetPprParameter(pulse.pprParamList,'t_pol',100)*1e-3;
            polarisationField = GetPprParameter(pulse.pprParamList,'b_pol',200)*1e-3;
            views_per_seg = GetPprParameter(pulse.pprParamList,'views_per_seg',4);
            acquisitionField = GetFrequencyInHz(pulse.pprParamList)/gammaH;
            npts = GetPprParameter(pulse.pprParamList,'no_samples',64);
            period = GetPeriodInSecond(pulse.pprParamList);
            trdd = GetPprParameter(pulse.pprParamList,'trdd',70)*1e-6;
            postpulseDelay = GetPprParameter(pulse.pprParamList,'taftfc2',0)*1e-3;
            echoTime = GetPprParameter(pulse.pprParamList,'te',0)*1e-3;
            effechoTime = GetPprParameter(pulse.pprParamList,'te_2',0)*1e-3;
            echoNumber = GetPprParameter(pulse.pprParamList,'no_echoes',1);
            lineNumber = GetPprParameter(pulse.pprParamList,'no_views',1);
            tramp = GetPprParameter(pulse.pprParamList,'tramp',1)*1e-6;
            tref = GetPprParameter(pulse.pprParamList,'tref',1)*1e-6;
            tdp = GetPprParameter(pulse.pprParamList,'tdp',1)*1e-6;
            acquisitionTime = (echoTime*views_per_seg + (echoTime*0.5) + (2*tramp) + (4e-6/2) + 2.2e-6 + 2e-6).*1.1; % add 10% time for good measure
            % recycle delay parameters
            recycleDelay = acquisitionTime*1.5+polarisationTime;
            % make the waveforms
            nView2 = GetPprParameter(pulse.pprParamList,'no_views_2',1);
            sliceNumber = GetPprParameter(pulse.pprParamList,'no_slices',1);
            no_experiments = GetPprParameter(pulse.pprParamList,'no_experiments',1);
            inversionPulseDuration = 4*1e-3;
           
            evolutionSwitch = GetPprParameter(pulse.pprParamList,'evol_on',1)&&...
                              ~isempty(pulse.waveformParam.evolutionField)&&...
                              GetPprParameter(pulse.pprParamList,'inv_obs_mod_level',1); 
            
            nView = GetPprParameter(pulse.pprParamList,'no_views',1);
            nView2 = GetPprParameter(pulse.pprParamList,'no_views2',1);
            iterationNumber = nView2*nView*sliceNumber/pulse.slicesPerWaveform;  
            if evolutionSwitch                 
                % list of evolution fields and times
                for index = 1:length(pulse.waveformParam.evolutionField)
                    evolutionTime = pulse.waveformParam.evolutionTime(index);
                    evolutionField = pulse.waveformParam.evolutionField(index);                    
                    if recycleDelay < (sum([polarisationTime,   inversionPulseDuration,   evolutionTime,      acquisitionTime])*1.2)
                        recycleDelay = sum([polarisationTime,   inversionPulseDuration,   evolutionTime,      acquisitionTime])*1.2;
                    end
                    pulse.waveformProfile = WaveformProfile(pulse.waveformProfile,...
                                                            [polarisationTime,   inversionPulseDuration,   evolutionTime,      acquisitionTime,   recycleDelay     ],...
                                                            [polarisationField,  acquisitionField,         evolutionField,     acquisitionField,  0                ],...
                                                            [0                   1                         0                   1                  0                ],...
                                                             pulse.dBdt,iterationNumber);
                end
            else
                pulse.waveformProfile = WaveformProfile(pulse.waveformProfile,...
                                                        [polarisationTime,   acquisitionTime,   recycleDelay     ],...
                                                        [polarisationField,  acquisitionField,  0                ],...
                                                        [0                   1                  0                ],...
                                                        pulse.dBdt,iterationNumber);
            end      
            % finally, sort out the number of averages
            pulse.waveformProfile.averageNumber = GetPprParameter(pulse.pprParamList,'no_averages',1);
            % test the waveform
            if ~TestSequence(pulse.waveformProfile)
                disp('Error: pulse sequence not safe.')
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
            end
            out = 1;
        end
        
%         function out = AcquisitionFunction(pulse,data)
%             % This sequence acquires two echoes per acquisition: the first
%             % one is an FID and can be used to calibrate the field, the
%             % second one is the image data. Each echo triggers this
%             % function once.  
%             
%             sigma = std(abs(pulse.scan.noiseSample));
%             rawData = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
%             echoNumber = data.Header{5}+1;
%             lineNumber = data.Header{2}+1;
%             expNumber = data.Header{6}+1;
%             nAcq = data.Header{9};  % number of acquisitions performed since the start of the scan
%             Nsample = data.Header{1};
% 
%             switch pulse.scan.GetEvoStatus
%                 case 'Setup'
%                 case 'Acquiring'
%                     fieldInd = pulse.waveformParam.indexList(expNumber,1);
%                     timeInd = pulse.waveformParam.indexList(expNumber,2);
%                     scanData = pulse.data{end}{2};
%                     scanData(:,lineNumber,echoNumber,timeInd,fieldInd) = rawData;
%                     pulse.data{end} = {pulse.pprParamList,scanData};
%                 otherwise
%             end
%             out = 1;
%         end
        
        function out = FinishAcquisition(pulse)
            if isequal(pulse.scan.acquisitionMode,'Acquiring')
                A = pulse.data;
                no_experiments = GetPprParameter(pulse.pprParamList,'no_experiments',1);
                sliceNumber = GetPprParameter(pulse.pprParamList,'no_slices',1);
                echoTime = GetPprParameter(pulse.pprParamList,'te',0)*1e-3;
                effechoTime = GetPprParameter(pulse.pprParamList,'te_2',0)*1e-3;
                views_per_seg = GetPprParameter(pulse.pprParamList,'views_per_seg',4);
                L = size(A,1);
                P = size(A,2);
                E = size(A,2)/no_experiments;
                S = E/sliceNumber;
                AA = zeros(L,L,sliceNumber,E);
                for e=1:no_experiments
                    for s=1:sliceNumber
                        AA(:,:,s,e) = A(:,1+(s-1)*S+(e-1)*E:S+(s-1)*S+(e-1)*E);
                    end %for slices
                end %for experiments
                %first the full kspaces

                teff = effechoTime./echoTime;
                tablename = [num2str(teff) '#' num2str(views_per_seg) '#' num2str(L) '.txt'];
                currentFolder = pwd;
                cd 'R:\CLSM\School of Medicine and Dentistry\Teams\Field-Cycling MRI\Structured_storage\Processing\Image Processing\PEtables'                                               %does the Pe reordering table exist?
                if exist(tablename)==0;
                    command = ['tgen ' tablename ' ' num2str(teff) ' ' num2str(views_per_seg)  ' ' num2str(L)]; %if not create it using the SMIS tool
                    dos(command);
                end
                petable = importdata(tablename);
                cd(currentFolder)
                petable = petable+1;


                dummy = zeros(L,L,sliceNumber,no_experiments);
                for e =1:2
                    for s=1:sliceNumber
                        for n=1:L
                            step = petable(n);
                            dummy(:,step,s,e)=AA(:,n,s,e);
                        end
                    end
                end
                data = dummy;
                im = fftshift(ifft2(data(:,:,1,1)));

%                 am_ref = mean(abs(navi(:,1)));
                
                figure('WindowStyle','Docked','Visible','on')
                imshow(abs(im)',[0 max(abs(im(:)))])
                disp(['Standard deviation: ' num2str(std(im(:)))])
            end
            out = 1;
        end
        
    end
    
end

