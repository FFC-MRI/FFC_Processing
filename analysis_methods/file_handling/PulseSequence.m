classdef PulseSequence < handle
    % PulseSequence is a class that helps with manipulation of pulse
    % sequence data and acquisition.
    
    properties
        pulseSequenceName           % name of the pulse sequence as it appears in the control window.
        dataFile                    % path and generic name of the data file.
        dataFileID                  % fid linking to the file (from fopen)
        pprFile                     % path and name of the ppr file 
        pprParamList                % list of the PPR parameters for quick access during in-scan corrections
        mrdFile                     % name of the MRD file saved by the EVO software
        data                        % raw data obtained from the scan
        dataProcessed               % data after post-processing
%         waveformParam               % parameters used for the waveform that cannot be handled correction by the PPR file (such as very low evolution fields)
        waveformProfile             % waveform object describing the field variation. 
        dBdt                        % time-variation of the field during the slopes (in T/s)
        useMode                     % 1 if setup-mode only, 2 if running-mode only, 3 if both can be used.
        defaultMode                 % default acquisition mode. 1 for setup, 2 for acquiring.
        autoProcess                 % if set to 1, start the image processing directly after acquisition.
        steppedAcquisition          % 1 if 'AcquisitionFunction' must be interrupted after each FID, 0 if it does not have to.
        autoStart                   % 1 if the sequence runs automatically at load time, 0 if the user has to start the sequence manually.
        slicesPerWaveform           % number of slices acquired for each waveform played. This is used to set the number of iterations automatically.
        preAquisitionDelay          % delay required by the EVO console, which can be ovelapped with the field stabilisation delay
        userData                    % additional data for the user, such as explanative text and so on.
    end
    
    properties (Access = private)
    end
    
    properties (Transient)
        scan                        % refence to the scan object that controls the scanner.
        editorHandle                % handle to the PPR editor activex control
    end
    
    methods
        
        function pulse = PulseSequence
            pulse.pulseSequenceName = '';
            pulse.dataFile = '';
            pulse.dataFileID = -1;
            pulse.pprFile = [];
            pulse.mrdFile = '';
            pulse.data = []; % contains a list of {PPR file, data from scan}
            pulse.dataProcessed = [];
            pulse.dBdt = 7;
            pulse.waveformProfile = WaveformProfile;
            pulse.waveformProfile.pulse = pulse;
            pulse.useMode = 3;
            pulse.defaultMode = 2;
            pulse.steppedAcquisition = 0;
            pulse.editorHandle = [];
            pulse.autoStart = 0;
            pulse.slicesPerWaveform = 1;
            pulse.pprParamList = {};
            pulse.preAquisitionDelay = 0;
            pulse.autoProcess = 1;
            pulse.userData.sequenceInfo = 'Write scan details here.';
            % you can set the default wavforms using the genrator object:
%             pulse.waveformProfile.generator.Tevo = 0; % this would remove the evolution period
        end
        
%% Handling the MRD file        
        
        function out = SetDataFile(pulse,fileName)
            % test that the file is a correct choice
            pathstr = fileparts(fileName);
            oldDir = cd;
            cd(pathstr)
            % set the file
            pulse.dataFile = fileName;
            cd(oldDir)
            out = 1;
        end
                
        function out = OpenDataFile(pulse)
            out = fopen(pulse.dataFile);
        end
        
        function out = TestDataFile(pulse)
            out = pulse.dataFileID;
        end
        
        function out = CloseDataFile(pulse)
            out = fclose(pulse.dataFileID);
        end
        
        
%% Handling the PPR file
        
        % list the PPR parameters for quick access during the scan
        function out = UpdatePprParameterList(pulse)
            if isempty(pulse.pprParamList)
                pulse.pprParamList = GetPprData(pulse.pprFile);
            elseif isobject(pulse.scan)
                if isequal(pulse,pulse.scan.pulse)
                    pulse.pprParamList = GetPprData(pulse.scan.evoHandle); % pull the EVO parameters if the sequence is loaded
                else
                    pulse.pprParamList = GetPprData(pulse.pprFile); % otherwise pull the PPR data file
                end
            else
                pulse.pprParamList = GetPprData(pulse.pprFile); % otherwise pull the PPR data file
            end
            
%             out = UploadPprFile(pulse);
            out = 1;
        end
        
        % uploads the ppr parameters to the EVO or to the PPR file
        function out = UploadPprFile(pulse)
            out = UploadPprFile(pulse.scan.evoHandle,pulse.pprParamList);
        end
        
%% Handling the pulse sequence data
        
        function out = SaveData(pulse)
            saveName = [pulse.dataFile(1:(end-3)) 'mat']; % replace '.MRD' with '.mat'
            save(saveName,pulse);
            out = 1;
        end
        
        % copy the pulse sequence into a new object. Useful to save data
        % after a scan if the sequence is to be used again
        function pulseCopy = CopySequence(pulse)
            pulseCopy = feval(class(pulse)); % create a new object with the same class
            fName = fields(pulseCopy); 
            fName(strcmp(fName,'scan')) = []; % remove the fields that should not be saved
            fName(strcmp(fName,'editorHandle')) = [];
            for ind = 1:length(fName)
                pulseCopy = setfield(pulseCopy,fName{ind},getfield(pulse,fName{ind})); %#ok<GFLD,SFLD>
            end
            pulseCopy.waveformProfile = Copy(pulse.waveformProfile); % waveformProfile is a handle, it needs particular care.
            pulseCopy.waveformProfile.pulse = pulseCopy;
        end
        
        % update the pulse object with a given list of evolution parameters
        function out = SetEvolutionParameter(pulse,fieldList,timeList)
            try
                validateattributes(fieldList,{'numeric'},{'>=',0,'<=',0.2});
                validateattributes(timeList,{'numeric'},{'>=',0,'<',2});
            catch
                out = 0;
                return
            end
            pulse.waveformProfile.generator.T1est = timeList;
            pulse.waveformProfile.generator.Bevo = fieldList;
            out = 1;
        end
        
        % determine the list of evolution fields and times manually
        function out = SelectEvolutionParameter(pulse)
            ManualEvolutionSelectionGui(pulse);
            out = 1;
        end
        
        % get the list of evolution fields and time using the organ
        % database
        function out = LookUpEvolutionParameter(pulse)
            out = 1;
        end
        
        % this function generates the waveform profile. It may be
        % overriden if needs be, though the default function is quite
        % general and should work in most cases.
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
                        inversionPulseDuration = 8000e-6;
                                            
                end
            end
            gradRampTime =  GetPprParameter(pulse.pprParamList,'tramp',1000)*1e-6*1.1;  % adding 10 % time for good measure
            inversionPulseDuration = inversionPulseDuration + gradRampTime*2;
            % The prepulse delay is already taken into accound within the
            % WaveformProfile object.
            % readout parameters
            acquisitionField = GetFrequencyInHz(pulse.pprParamList)/gammaH;
            npts = GetPprParameter(pulse.pprParamList,'no_samples',64);
            period = GetPeriodInSecond(pulse.pprParamList);
            trdd = GetPprParameter(pulse.pprParamList,'trdd',70)*1e-6;
            postpulseDelay = GetPprParameter(pulse.pprParamList,'taftfc2',0)*1e-3;
            echoTime = GetPprParameter(pulse.pprParamList,'te',0)*1e-3;
            echoNumber = GetPprParameter(pulse.pprParamList,'no_echoes',1);
            acquisitionTime = (gradRampTime*2 + npts*period + trdd + postpulseDelay + echoTime*echoNumber)*1.1; % add 10% time for good measure
            sliceNumber = GetPprParameter(pulse.pprParamList,'no_slices',1);
            nView = GetPprParameter(pulse.pprParamList,'no_views',1);
            nView2 = GetPprParameter(pulse.pprParamList,'no_views2',1);
            iterationNumber = nView2*nView*sliceNumber/pulse.slicesPerWaveform;  
            % make the waveforms
            pulse.waveformProfile.generator.Bread = acquisitionField;
            pulse.waveformProfile.generator.Tread = acquisitionTime;
            pulse.waveformProfile.generator.inversionTime = inversionPulseDuration;
            pulse.waveformProfile.generator.iterationNumber = iterationNumber;
            pulse.waveformProfile.generator.averageNumber = GetPprParameter(pulse.pprParamList,'no_averages',1);
            pulse.waveformProfile.generator.waveform.pulse = pulse;
                 
            [okflag,pulse.waveformProfile.generator] = Update(pulse.waveformProfile.generator);
            if ~okflag
                error('Invalid waveform. Cannot proceed with the sequence.')
            end
        end
        
        function out = UploadWaveform(pulse)
            out =  UploadWaveform(pulse.waveformProfile);
        end
        
        % this function is used before the scan to allow synchronising the
        % sequence with some of the scan parameters, such as the coil. It
        % should be overriden by derived classes, if necessary.
        function out = SyncPulseSequence(pulse) %#ok<*MANU,*INUSD>
%             pulse.mrdFile = pulse.scan.evoHandle.DataFile;
            out = 1;
        end
        
        
        
%% Functions called for the scan (while the pulse sequence is loaded)

        % function called just before starting an acquisition
        function out = PreparePulseSequence(pulse)
            pulse.data = zeros(1,1,1,1,1,1,1,1,1);    % prepare a slot for the data acquisition
            switch pulse.scan.acquisitionMode
                case 'Acquiring'  % acquiring
                    
                case 'Setup'
            end
            out = 1;
        end
        
        % Function that starts the pulse sequence. 'mode' can be 'Setup' or
        % 'Acquiring'
        function out = RunPulseSequence(pulse)
            if pulse.steppedAcquisition
                if StartSteppedAcquisition(pulse.scan)
                    while NextStep(pulse.scan) % continues stepping until the EVO console stops
                        while (isequal(ScanStatus(pulse.scan),'Stepping')&&IsEvoAcquiring(pulse.scan))   % wait for the step to complete or for the user to stop the setup
                            pause(0.01)
                        end
                    end
                    out = StopSteppedAcquisition(pulse.scan);
                else
                    disp('Missed step in the acquisiton.')
                    out = 0;
                end
            else
                out = ContinuousAcquisition(pulse.scan);
            end
        end
        
        % method called during the scan when the ScanData event is
        % triggered.
        % Classes derived from this class
        % should override this function to include a processing function
        % inside the acquisition loop.
        % Beware that this function is designed to act on a ScanControl
        % object.
        % Be careful when using this function: you cannot change the PPR
        % parameters on the EVO system while it is running! so do not use
        % SetPprData or the likes
        function out = AcquisitionFunction(pulse,data)
            % Use this in a subclass to run this code:
%             out = AcquisitionFunction@PulseSequence(pulse,data);

            % The code below sorts out the data correctly accoding to the
            % various loops and receivers.
            acqLength = data.Header{1};
            lineNumber = data.Header{2}+1;
            views2Number = data.Header{3}+1;
            sliceNumber = data.Header{4}+1;
            echoNumber = data.Header{5}+1;
            expNumber = data.Header{6}+1;
            averageNumber = data.Header{7};
            channelNumber = data.Header{8}+1;
            nAcq = data.Header{9};  % number of acquisitions performed since the start of the scan
            % reinitialise the data if it is the very first acquisition
            if nAcq==1
                pulse.data = zeros(acqLength,1,1,1,1,1,1,1,1);
            end
            % store the data received
%             if isequal(pulse.scan.acquisitionMode,'Acquiring') || pulse.steppedAcquisition
            try
                fieldInd = pulse.waveformProfile.waveList{expNumber}.indexBevo;
                timeInd = pulse.waveformProfile.waveList{expNumber}.indexTevo;
            catch
                fieldInd = 1;
                timeInd = expNumber;
            end
            pulse.data(1:size(data.DataBuffer,2),...
                       lineNumber,...
                       views2Number,...
                       sliceNumber,...
                       echoNumber,...
                       timeInd,...
                       fieldInd,...
                       averageNumber,...
                       channelNumber) = double(data.DataBuffer(1,:)) + 1i*double(data.DataBuffer(2,:));
%             end
            % make a pause if the acquisition is stepped (stop on the first
            % dataset received to make this responsive)
            if pulse.steppedAcquisition && channelNumber==1 && echoNumber==1
                pulse.scan.scannerStatus = 'Paused';
            end
%             try
%                 disp(pulse.scan.evoHandle.ProgressTask)
%                 disp(hex2dec(pulse.scan.evoHandle.Progress))
%             end
            out = 1;
        end
        
        % this function is a post-acquisition function that gets called
        % before the pulse sequence is unloaded. It does not use parallel
        % processing so it should only be used for quick processing such as
        % data display
        function out = FinishAcquisition(pulse)
            % stop the scanner if it is in stepping mode
            out = 1;
            if pulse.steppedAcquisition
                out = AbortEvo(pulse.scan);
            end
        end
        
%% Post-processing section: functions called after the pulse is unloaded from the scan
        
        % the next function processes the raw data file and generates the
        % data. This should be overriden by the classes derived.
        % Beware that this function may be proecssed in the parallel pool,
        % in which case the handle 'pulse' will NOT be updated after the
        % function is finished, unless it is passed back to the calling
        % function like a normal, non-handle variable.
        % (See the ControlScan function 'LaunchParallelProcessing')
        function [out,pulse] = ProcessData(pulse)
            pulse.dataProcessed = pulse.data; % by default, no operation is done
            out = 1;
        end
        
        % method to show the data in a figure, launched automatically after
        % data processing
        function out = ShowData(pulse)
%             plot(pulse.data{end}{2})
%             plot(pulse.data)
            out = 1;
        end
        
    end
    
end