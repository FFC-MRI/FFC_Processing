classdef WaveformGenerator
    %WaveformGenerator Class used to generate waveform objects from
    %elementary data such as evolution field, evolution time, inversion
    %pulse, half-Fourier space etc...
    %   The WaveformProfile objects contain an instance of this class that
    %   is used to generate the waveform profiles just before the scan
    %   begins. You may adapt the generator object to the pulse sequence
    %   and use particular instances by calling them within the
    %   PulseSequence object when you create the pulse object.
    %   Each value in the properties may be a single number or an array. If
    %   it is an array then its value is changed for each iteration.
    
    properties
        name = 'default';   % name of the generator (should be changed for each type of waveform generator)
        Bpol = 0.2;     % polarisation field, in T
        Tpol = 0.2;     % polarisation time, in sec
        Bevo = 0.1;     % evolution field, in T
        Tevo = 0;       % evolution time in sec. This may be empty or set to 0 if 
                        % no evolution is used, an array of only one
                        % evolution is required (two-point method, for
                        % instance) or a matrix if several evolutions are
                        % required for each field. Lines: fields, Cols:
                        % times
        iterationNumber = 1;        % number of times each waveform is repeated (number of lines in the image if all k-space lines are acquired)
        averageNumber = 1;          % number of repetitions for each profile
        Bread = 0.198;          % readout field, in T
        Bfid = 0;               % transmit field for double resonance experiment acquisition
        Tread = 0.05;           % readout time, in sec
        readDelay = 0.015;      % delay before readout, added to the system delay for field stabilisation.
        Tcool = 0.1;            % cooldown time. If lower than the minimum required, this value is set to the minimum allowed.
        Tinv = 0;               % duration of the inversion pulse, listed here for reference
        
        doubleResonance = false;    % flag to use double resonance irradiation during evolution.
        waveform = [];              % WaveformProfile object that called this generator. This reference is used for convenience.
        waveformBlank = [];         % empty waveform used for k-space filling, for partial acquisitions
        prePulseDelay = 0;          % additional delay before the RF pulses, on top of the one set by the waveform object. Can be negative.
        % additional parameters to make it easier to build the lists of
        % Bevo and Tevo:
        BevoStart = 0.2;          % starting value for the Bevo list
        BevoEnd = 0.2;      % end value for the Bevo list
        NBevo = 1;              % number of evolution fields
        T1est = 0;           % estimation of T1 for each field
        TevoStartFactor = 2.5;      % multiplication factor to generate the list of Tevo values
        TevoEndFactor = 0.2;
        NTevo = 1;              % number of evolution times per field
        timeAxis = 'log';       % accquisition mode for the evolution time, may be 'lin' or 'log'
        qpSampling = 'none';    % sampling quality of the QP
        organReference = [];    % reference object for the tissue investigated
        indexList = [];         % list of correspondance indexes for Bevo and Tevo, for convenience
        undersamplingMethod = 'none';       % undersampling method to use for the k-space filling
        undersamplesize = 1;
        accelerationfactor =1;
        kspacelineaquired = [];
        roundedCorner = 1;      % rounds the edges of the field profile
        inversionTime = 0;      % duration of the inversion pulse in the case of inversion recovery
        inversionFlag = false;  % array of flags for inversion recovery experiments
        
        % parameters for non-adiabatic transition experiments
        naVector = {[1 0 0]};     % orientation of the NA vector
        naTriggerDown = false;               % boolean, set the controls for non-adiabatic transitions for the down slope.
        naTriggerUp = false;               % boolean, set the controls for non-adiabatic transitions for the up slope.
        naDelayDown = 0;            % delay in sec between the end of the down ramp and the desactivation of the NA coil.
        naDurationDown = 0;         % duratio of the activation of the NA coil during the down ramp
        naDelayUp = 0;              % delay in sec between the start of the up ramp and the activation of the NA coil.
        naDurationUp = 0;         % duratio of the activation of the NA coil during the up ramp  
        waveformIndex = 0;          % selection index for the waveforms
    end
    
    properties (Access = protected)
        dBdtMax = 15;               % max field variation (T/s)
        dBdtMin = 0.01;
        dBdtComfort = 4;            % field slope used to start and stop the profiles, to minimise patient discomfort (T/s)
        dBdtPerf = 10;              % field slope used  to minimise magnetisation losses during the experiment (T/s)
        delatchAmplitude = 3.6e-3;  % amplitude of the delatch element (T)
        delatchDuration = 1e-3;     % in s
        delatchDelay = 1e-3;        % delay between the delatch and the main waveform
        prepulseDelay = 15e-3;      % delay to stabilise the field before the RF pulse is sent
        maxField = 0.2;             % maximum field value (T)
        minField = 0;               % minimum field value (T)
        IECOMaxPowerOut             % to be confirmed, this can be used to test if the sequence risks overheating the amps
        minDeltaT = 1e-5;           % minimum duration used to flag a TTL line without delay (to avoid interpolation problems with the cDAQ)
        currentPerTesla = 1800/0.2; % in A/T (needs fine tuning)
        voltPerTesla = 100/0.2;     % in V/T (needs confirmation)
        coilResistance = 83e-3/3;   % overall resistance of the magnet
        maxAveragePower = 60e3;     % in W, maximum power dissipation of the chiller unit. Initially 60e3
        cDaqDelay = 0.35;           % additional cooldown delay due to the cDAQ program
    end
    
    methods
        function wavegen = WaveformGenerator(varargin)
            % varargin is a waveform profile
            % this constructor method adds waveforms to the object. If the
            % object does not exist, it creates one.
            if nargin == 1
                wavegen.waveform = varargin{1};
            else
                wavegen.waveform = WaveformProfile(wavegen);
            end
            wavegen.waveformBlank = WaveformProfile(wavegen);
        end
        
        function out = ShowPrivateProperty(wavegen,paramName)
            try 
                out = getfield(wavegen,paramName); %#ok<GFLD>
            catch ME
                error('Field does not exist');
            end
        end
        
        % transfer the value of the parameters from the old to the new
        % waveform generator
        function wavegen = TransferParameter(wavegen,oldwavegen)
            proplist = properties(oldwavegen);
            for ind = 1:length(proplist)
                if isprop(wavegen,proplist{ind})
                    if ~isequal(proplist{ind},'name')
                        wavegen = setfield(wavegen,proplist{ind},getfield(oldwavegen,proplist{ind})); %#ok<GFLD,SFLD>
                    end
                end
            end
        end
        
        % make a copy of the waveform
        function waveCopy = Copy(wavegen)
            waveCopy = feval(class(wavegen)); % create a new object with the same class
            fName = fields(waveCopy);
            for ind = 1:length(fName)
                waveCopy = setfield(waveCopy,fName{ind},getfield(wavegen,fName{ind})); %#ok<GFLD,SFLD>
            end
        end
        
        % set the new generator to a waveform profile object
        function wavegen = TransferGenerator(wavegen,oldwavegen)
            wavegen = TransferParameter(wavegen,oldwavegen);
            if ~isempty(oldwavegen.waveform)
                oldwavegen.waveform.generator = wavegen;
            end
        end
        
        % provide an estimation of the duration of the entire experiment
        function t = EstimateDuration(wavegen)
            t = 0;
            for i = 1:length(wavegen.waveform.waveList)
                t = t+wavegen.waveform.waveList{i}.time(end)*wavegen.waveform.waveList{i}.iterationNumber;
            end
            t = t*wavegen.waveform.averageNumber;  % time in seconds
        end
        
        % add some points over the quadrupolar peaks
        function wavegen = ProbeQP(wavegen)
            switch lower(wavegen.qpSampling)
                case 'none'
                case '1-pt'
                    wavegen.Bevo = [wavegen.Bevo; 0.0658];
                case 'very low'
                    wavegen.Bevo = [wavegen.Bevo; 0.0939; 0.0658; 0.0352];
                case 'low'
                    wavegen.Bevo = [wavegen.Bevo; 0.0939; 0.0658; 0.0564; 0.0470; 0.0352];
                case 'medium'
                    wavegen.Bevo = [wavegen.Bevo; logspace(log10(3.5),log10(1.6),10)'*1e6/gammaH];
                case 'high'
                    wavegen.Bevo = [wavegen.Bevo; logspace(log10(3.5),log10(1.6),15)'*1e6/gammaH];
                case 'very high'
                    wavegen.Bevo = [wavegen.Bevo; logspace(log10(3.5),log10(1.6),30)'*1e6/gammaH];
            end
        end
        
        function wavegen = UpdateField(wavegen)
            % make the list of evolution fields
            wavegen.Bevo = logspace(log10(wavegen.BevoStart),log10(wavegen.BevoEnd),wavegen.NBevo)'; % the first dimension is set to the number of fields
            wavegen = ProbeQP(wavegen); % adds the sampling for the QP
            wavegen.Bevo = sort(wavegen.Bevo,'descend');
            if length(wavegen.Bevo)~=length(wavegen.inversionFlag)
                wavegen.inversionFlag = false(size(wavegen.Bevo)); % re-generate the inversion flag list   
            end
        end
        
        % Provide T1 estimation at each field from the tissue reference
        % object
        function wavegen = UpdateT1Estimation(wavegen)
            if ~isempty(wavegen.organReference)
                try
                    wavegen.T1est = 1./wavegen.organReference.dispersionModel(wavegen.Bevo);
                catch
                    error('Invalid tissue model')
                end
            else
                wavegen.T1est = wavegen.T1est(1)*ones(size(wavegen.Bevo));
            end
            wavegen.T1est = round(wavegen.T1est*1e6)/1e6; % limits the precision to 1us
        end
        
        % update the list of evolution times
        function wavegen = UpdateTime(wavegen)
            wavegen.Tevo = zeros(length(wavegen.Bevo),wavegen.NTevo);
            for indField = 1:size(wavegen.Bevo,1)
                if wavegen.NTevo == 1
                    wavegen.Tevo(indField,:) = wavegen.T1est(indField);
                else
                    switch wavegen.timeAxis
                        case {'lin','Lin','LIN'}
                            wavegen.Tevo(indField,:) = linspace(wavegen.T1est(indField)*wavegen.TevoStartFactor,...
                                                                wavegen.T1est(indField)*wavegen.TevoEndFactor,...
                                                                wavegen.NTevo);
                        case {'log','Log','LOG'}
                            wavegen.Tevo(indField,:) = logspace(log10(wavegen.T1est(indField)*wavegen.TevoStartFactor),...
                                                                log10(wavegen.T1est(indField)*wavegen.TevoEndFactor),...
                                                                wavegen.NTevo);
                    end
                end
            end
            wavegen.Tevo = round(wavegen.Tevo*1e6)/1e6; % limits the precision to 1us
        end
        
        % This function updates the various waveform parameters according
        % to NBevo, T1est and so on. It allows desiging many waveforms
        % more easily.
        function [out,wavegen] = Update(wavegen,estimationFlag)
            out = UpdatePprParameterList(wavegen.waveform.pulse);
            if nargin < 2
                estimationFlag = true; % give a possibility to bypass the automatic estimation
            end
            wavegen = UpdateField(wavegen);
            if estimationFlag
                wavegen = UpdateT1Estimation(wavegen);
            end
            wavegen = UpdateTime(wavegen);
            % generate the waveforms and populate the waveform object
            wavegen = MakeWaveform(wavegen);
            wavegen = UnderSample(wavegen);
            wavegen = UpdatePulse(wavegen);
            out = TestSequence(wavegen);
        end
                        
        % generate the waveforms
        function wavegen = MakeWaveform(wavegen)
            % This is the basic waveform used by default for waveform
            % generators. It is a simple polarisation - evolution
            % - detection pattern using the parameters provided. 
                                    
            % reinitialise the fields
            wavegen.waveform.waveList = {};
            wavegen.waveformBlank.waveList = {};
            wavegen.indexList = [];
            wavegen.waveformIndex = 0;
                
            % generate the waveforms
            for indField = 1:length(wavegen.Bevo)
                for indTevo = 1:size(wavegen.Tevo,2)
                    % select the next index and create the list
                    wavegen.waveformIndex = wavegen.waveformIndex +1;
                    ind = wavegen.waveformIndex;
                    % store the index of the evolution field and time into
                    % an array:
                    wavegen.waveform.waveList{ind}.indexBevo = indField;
                    wavegen.waveform.waveList{ind}.indexTevo = indTevo;
                    wavegen.indexList(ind,:) = [indField indTevo];
                    
                    % Start the waveform shape here:
                    wavegen = InitialiseWaveform(wavegen);                    
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Bevo = wavegen.Bevo(indField);
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Tevo = wavegen.Tevo(indField,indTevo);
                    
                    % make the waveform segment by segment
                    wavegen = Delatch(wavegen);   % adds the delatch signal for the IECO amps
                    % polarisation period:
                    if wavegen.Tpol~=0
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bpol);  % go to polarisation 
                        wavegen = Plateau(wavegen,wavegen.Tpol);      % stay at the polarisation field for the polarisation time
                    end
                    % now going to evolution field:
                    if wavegen.Tevo(indField,indTevo) ~= 0
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bevo(indField));
                        wavegen = Plateau(wavegen,wavegen.Tevo(indField,indTevo));
                    end
                    % and finally go to readout                    
                    wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bread);
                    wavegen = PlateauWithPulse(wavegen,wavegen.Tread,wavegen.readDelay,0);
                    % ends with the cooling time
                    wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
                    if wavegen.Tcool < minCoolTime(wavegen)
                        wavegen = Plateau(wavegen,minCoolTime(wavegen)*1.05);
                    else
                        wavegen = Plateau(wavegen,wavegen.Tcool);
                    end
%                     wavegen.waveform
                end
            end
            % find out the number of experiments to program the EVO. This
            % has to be done before the undersampling.
            wavegen.waveform.experimentNumber = wavegen.waveformIndex;
            wavegen.waveform.averageNumber = wavegen.averageNumber;
            wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber = wavegen.iterationNumber;
%             wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber = GetPprParameter(wavegen.waveform.pulse.pprParamList,'no_views',1);
        end        
        
        % undersample the k-space. kSpaceLine is an array of booleans that
        % is set to 'true' if the line is measured, 'false' if it is
        % skipped
        function wavegen = UnderSample(wavegen)
            
            % Generates a list of booleans that select the k-space lines to
            % remove
            no_samples = GetPprParameter(wavegen.waveform.pulse.pprParamList,'no_samples',1);
            switch wavegen.undersamplingMethod
                case 'none'
                    wavegen.kspacelineaquired = 1:no_samples;
                    return
                case 'Partial Fourier'
                    wavegen.kspacelineaquired = UndersampleMask(no_samples,'PF');
                case 'Compressed Sensing'
                    wavegen.kspacelineaquired = UndersampleMask(no_samples,'CS');
            end
            % new way to do: only take 75% of the phase encoding lines
            nlines = sum(wavegen.kspacelineaquired);
            wavegen.waveform.pulse.pprParamList = SetPprData(wavegen.waveform.pulse.pprParamList,'effective_no_views',no_samples);
            wavegen.waveform.pulse.pprParamList = SetPprData(wavegen.waveform.pulse.pprParamList,'no_views',nlines);
            for nexp = 1:length(wavegen.waveform.waveList)
                wavegen.waveform.waveList{nexp}.iterationNumber = nlines;
            end
%             % runs through the list of k-space line to acquire and select
%             % the correct waveform with corresponding number of iterations:
%             completeWaveList = {};
%             for indexp = 1:wavegen.waveform.experimentNumber % unwrap each experiment individually (not optimal, but does not matter)
%                 lineIndex = 1;
%                 while lineIndex <= length(kSpaceLine)
%                     % groups all the lines with the same acquisition protocol
%                     % as the current one to gain time on data transmission
%                     bundle = lineIndex;
%                     for ind = (lineIndex+1):length(kSpaceLine)
%                         if kSpaceLine(lineIndex) == kSpaceLine(ind)
%                             bundle = [bundle ind]; %#ok<AGROW>
%                         else
%                             break
%                         end
%                     end
%                     % list the protocol corresponding to the current bundle of
%                     % k-space lines
%                     if kSpaceLine(lineIndex) % case when we make the full acquisition
%                         completeWaveList = [completeWaveList wavegen.waveform.waveList(indexp)]; %#ok<AGROW>
%                     else % case of a blank acquisition for speed-up
%                         completeWaveList = [completeWaveList wavegen.waveformBlank.waveList(indexp)]; %#ok<AGROW>
%                     end
%                     completeWaveList{end}.iterationNumber = length(bundle);
%                     if bundle(end) < length(kSpaceLine)
%                         lineIndex = ind; % pass on to the next bundle
%                     else
%                         lineIndex = length(kSpaceLine)+1; % avoid the case when the last line is the same as the previous one, this causes problems
%                     end
%                 end
%             end
%             % finally, update the waveform object with the undersampled
%             % version:
%             wavegen.waveform.waveList = completeWaveList;  
        end
        
        % this function holds custom scripts that run after the generator
        % has finished creating the waveforms, it can be used to update
        % some fields in the PPL file.
        function wavegen = UpdatePulse(wavegen)
            % note: the pulse object can be accessed at
            % wavegen.waveform.pulse
            
            if ~isobject(wavegen.waveform.pulse)
                return
            end
                        
            % set the number of experiments in the pulse list to fit the
            % number of experiments to be done by the EVO.
            SetPprData(wavegen.waveform.pulse,'no_experiments',wavegen.waveform.experimentNumber);
            SetPprData(wavegen.waveform.pulse,'no_averages',wavegen.waveform.averageNumber);
            
            % copy the waveform parameters into the PPR file
            if ~isempty(wavegen.Tevo)
                times = zeros(1,50);
                times(1:length(wavegen.Tevo(:))) = wavegen.Tevo(:)';
                SetPprData(wavegen.waveform.pulse,'t_evol',  times(1:50)*1e3);
            end
            if ~isempty(wavegen.Bevo)
                fields = zeros(1,50);
                fields(1:length(wavegen.Bevo)) = wavegen.Bevo(:)';
                SetPprData(wavegen.waveform.pulse,'b_evol',fields(1:50)*1e3);
            end
            if ~isempty(wavegen.Bpol)
                val = wavegen.Bpol;
                SetPprData(wavegen.waveform.pulse,'b_pol',val*1e3);
            end
            if ~isempty(wavegen.Tpol)
                val = wavegen.Tpol;
                SetPprData(wavegen.waveform.pulse,'t_pol',val*1000);
            end
            if ~isempty(wavegen.Tcool)
                val = wavegen.Tcool;
                SetPprData(wavegen.waveform.pulse,'down_time',val*1000);
            end
            % by default, no inversion
            SetPprData(wavegen.waveform.pulse,'inv_obs_mod_level', 0);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic waveform shapes 
                
        function wavegen = InitialiseWaveform(wavegen)
            if wavegen.waveformIndex == 0
                wavegen.waveformIndex = 1;
            end
             % initialise the fields for the current waveform
            wavegen.waveform.waveList{wavegen.waveformIndex}.time = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.field = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber = wavegen.iterationNumber;
            wavegen.waveform.waveList{wavegen.waveformIndex}.naVector = [0 0 0];
            % do the same for the blank waveform
            wavegen.waveformBlank.waveList{wavegen.waveformIndex}.iterationNumber = wavegen.iterationNumber;
            wavegen.waveformBlank.waveList{wavegen.waveformIndex}.time = 0;
            wavegen.waveformBlank.waveList{wavegen.waveformIndex}.field = 0;
            wavegen.waveformBlank.waveList{wavegen.waveformIndex}.naTrigger = 0;
            wavegen.waveformBlank.waveList{wavegen.waveformIndex}.evoTrigger = 0;
            wavegen.waveformBlank.waveList{wavegen.waveformIndex}.naVector = [0 0 0];
        end
        
        function wavegen = Plateau(wavegen,duration)
            % waveform = Plateau(waveform,duration)
            % A basic field plateau, without acquisition. Nothing exciting.
            
            % test the inputs
            if ~isobject(wavegen.waveform)
                error('First input must be a WaveformProfile waveformect')
            end
            if duration < 0
                error('Negative duration.')
            end
            
            if duration>0
                % set the values
                wavegen.waveform.waveList{wavegen.waveformIndex}.time(end+1) = wavegen.waveform.waveList{wavegen.waveformIndex}.time(end) + duration;
                wavegen.waveform.waveList{wavegen.waveformIndex}.field(end+1) = wavegen.waveform.waveList{wavegen.waveformIndex}.field(end);
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(end+1) = 0;
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(end+1) = 0;
            end
        end
        
        function wavegen = PlateauWithPulse(wavegen,acqDuration,acquitionDelay,finishDelay)
            % waveform = PlateauWithPulse(waveform,acqDuration,acquitionDelay,finishDelay)
            % Maintains the field and sets the acquisition line after a
            % delay acquitionDelay. The field is maintained after the
            % acquisition is over for a duration finishDelay.
            
            % test the inputs
            if ~isobject(wavegen.waveform)
                error('First input must be a WaveformProfile waveformect')
            end
            if acqDuration < 0
                error('Negative duration.')
            end
            if acquitionDelay <= 0
                acquitionDelay = 1e-6;
                if acquitionDelay < 0
                    disp('Warning: negative pre-acquisition delay set to 0.')
                end
            end
            if finishDelay <= 0
                finishDelay = 1e-6;
                if finishDelay < 0
                    disp('Warning: negative post-acquisition delay set to 0.')
                end
            end
                
            if acqDuration>0
                % set the values for the acquisition waveform
                duration =                                                    [                                                            acquitionDelay       wavegen.minDeltaT    acqDuration          wavegen.minDeltaT   ]; 
                wavegen.waveform.waveList{wavegen.waveformIndex}.field =      [wavegen.waveform.waveList{wavegen.waveformIndex}.field       wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)  wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)  wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)   wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)];
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger =  [wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger   0                    0                    0                     0                  ];
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger = [wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  0                    1                    1                     0                  ];
                wavegen.waveform.waveList{wavegen.waveformIndex}.time =       [wavegen.waveform.waveList{wavegen.waveformIndex}.time  wavegen.waveform.waveList{wavegen.waveformIndex}.time(end)+cumsum(duration) ]; 
                % make the corresponding blank waveform
                wavegen.waveformBlank.waveList{wavegen.waveformIndex}.time =       [wavegen.waveformBlank.waveList{wavegen.waveformIndex}.time  wavegen.waveformBlank.waveList{wavegen.waveformIndex}.time(end)+cumsum(duration) ]; 
                wavegen.waveformBlank.waveList{wavegen.waveformIndex}.field =      zeros(size(wavegen.waveformBlank.waveList{wavegen.waveformIndex}.time));
                wavegen.waveformBlank.waveList{wavegen.waveformIndex}.naTrigger =  zeros(size(wavegen.waveformBlank.waveList{wavegen.waveformIndex}.time));
                wavegen.waveformBlank.waveList{wavegen.waveformIndex}.evoTrigger = [wavegen.waveformBlank.waveList{wavegen.waveformIndex}.evoTrigger    0  1  1  0];
                
            end
            if finishDelay > 0
                wavegen = Plateau(wavegen,finishDelay);
            end
        end
        
        
        function wavegen = RampConstantSlope(wavegen,slope,fieldEnd)
            % waveform = RampConstantSlope(waveform,slope,fieldEnd)
            % This method appends a slope to the list of elements, which is
            % defined by its slope (absolute value)
            
            % test the inputs
            if ~isobject(wavegen.waveform)
                error('First input must be a WaveformProfile object.')
            end
            if (slope > wavegen.dBdtMax)||(slope < wavegen.dBdtMin)
                error('The field time-derivate is out of boundaries.')
            end
            if (fieldEnd > wavegen.maxField)||(fieldEnd < wavegen.minField)
                error('The field is out of boundaries.')
            end
            
            % add the ramp to the current waveform object
            rampDuration = abs(wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)-fieldEnd)/slope;
            if rampDuration>0
                wavegen.waveform.waveList{wavegen.waveformIndex}.time(end+1) = wavegen.waveform.waveList{wavegen.waveformIndex}.time(end) + rampDuration;
                wavegen.waveform.waveList{wavegen.waveformIndex}.field(end+1) = fieldEnd;
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(end+1) = 0;
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(end+1) = 0;
            end
        end
        
        function wavegen = RampConstantTime(wavegen,rampDuration,fieldEnd)
            % waveform = RampConstantTime(waveform,rampDuration,fieldEnd)
            % This method appends a slope to the list of elements, which is
            % defined by its duration
            
            % test the inputs
            if ~isobject(wavegen.waveform)
                error('First input must be a WaveformProfile waveformect')
            end
            if (fieldEnd > wavegen.waveform.maxField)||(fieldEnd < wavegen.waveform.minField)
                error('The field is out of boundaries.')
            end
            dBdt = abs(wavegen.waveform.waveList{wavegen.waveformIndex}.field(end) - fieldEnd)/rampDuration;
            if (dBdt > wavegen.waveform.dBdtMax)
                error('The slope duration is too short.')
            end
            if (dBdt < wavegen.waveform.dBdtMin)
                error('The slope duration is too long.')
            end
            
            % add the ramp to the current waveformect
            if rampDuration>0
                wavegen.waveform.waveList{wavegen.waveformIndex}.time(end+1) = wavegen.waveform.waveList{wavegen.waveformIndex}.time(end) + rampDuration;
                wavegen.waveform.waveList{wavegen.waveformIndex}.field(end+1) = fieldEnd;
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(end+1) = 0;
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(end+1) = 0;
            end
        end
        
        function wavegen = StartWaveform(wavegen,delay,field)
            wavegen = Plateau(wavegen,delay);
            wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,field);
        end
        
        
        
        function wavegen = StopWaveform(wavegen,delay)
            wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
            wavegen = Plateau(wavegen,delay);
        end
        
        function wavegen = Delatch(wavegen)
            wavegen = StartWaveform(wavegen,0,wavegen.delatchAmplitude);
            wavegen = Plateau(wavegen,wavegen.delatchDuration);
            wavegen = StopWaveform(wavegen,wavegen.delatchDelay);
            wavegen = Plateau(wavegen,wavegen.delatchDuration);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Safety testing functions

        function energy = CalculateHeatEnergy(wavegen,ind)
            % this calculation only takes into account the power deposition
            % due to Joule heating. It processes the power deposition in
            % each of the three windings.
            if nargin == 1
                ind = wavegen.waveformIndex;
            end
            current = wavegen.waveform.waveList{ind}.field*wavegen.currentPerTesla;
            resistor = wavegen.coilResistance;
            time = wavegen.waveform.waveList{ind}.time;
            energy = trapz(time,current.^2)*resistor;
        end
        
        function IecoCoolTime = CalculateIecoCoolTime(wavegen,ind)
            % this functions makes sure that the IECO amplifiers have a 50%
            % duty ratio when operating at max power, while preserving DC
            % mode below 0.13 T.  
            if nargin == 1
                ind = wavegen.waveformIndex;
            end
            B0 = wavegen.waveform.waveList{ind}.field;
            time = wavegen.waveform.waveList{ind}.time;
            DR = @(b0) (-1/0.14*b0+0.27/0.14).*(b0>0.13) + (b0<=0.13);
            t = linspace(time(1),time(end),1000);
            step = t(2)-t(1);
            b = interp1(time,B0,t);
            IecoCoolTime =  sum(DR(b)*step);
        end
        
        function energy = CalculateCoolingEnergy(wavegen,ind)
            % derives the maximum energy evacuated by the chillers during
            % one cycle, supposing a steady state.
            if nargin == 1
                ind = wavegen.waveformIndex;
            end
            energy = (wavegen.waveform.waveList{ind}.time(end) + wavegen.cDaqDelay)*wavegen.maxAveragePower;
%             energy = energy/2;
        end
        
        % find the minimum cooling time for the currently selected waveform
        function minCoolTime = minCoolTime(wavegen)
            energyHeat = CalculateHeatEnergy(wavegen);
%             minOverallCoolTime = energyHeat/wavegen.maxAveragePower;
            MagnetCoolTime = energyHeat/wavegen.maxAveragePower;
%             minCoolTime = minOverallCoolTime...
%                           - wavegen.waveform.waveList{wavegen.waveformIndex}.time(end)...
%                           - wavegen.cDaqDelay; % min cooling time for the chiller
            IecoCoolTime = CalculateIecoCoolTime(wavegen);
            minCoolTime = max(MagnetCoolTime,IecoCoolTime);
            minCoolTime = minCoolTime - wavegen.cDaqDelay;
            if minCoolTime < 0
              disp(['Cooling time exceed the minimum allowed by ' num2str(abs(round(minCoolTime))) ' sec (which is good).'])
              minCoolTime = 0; % case when the cooling is already sufficient
            end
        end
        
        function out = TestSequence(wavegen)
            % out = TestSequence(waveform)
            % Test the sequence for problems, such as overheating,
            % correct length for all the fields, NaN, last value of the
            % field, etc...
            out = 1;            
            if isempty(wavegen.waveform.averageNumber)
                disp('Error: number of averages is not defined.')
                out = 0;
            elseif ~isnumeric(wavegen.waveform.averageNumber)
                disp('Error: number of averages is not numeric.')
                out = 0;
            elseif ~isfinite(wavegen.waveform.averageNumber)
                disp('Error: infinite number of averages.')
                out = 0;
            elseif wavegen.waveform.averageNumber<1
                disp('Error: negative number of averages.')
                out = 0;
            end
            if isempty(wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber)
                disp('Error: number of iteration is not defined.')
                out = 0;
            elseif ~isnumeric(wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber)
                disp('Error: number of iteration is not numeric.')
                out = 0;
            elseif ~isfinite(wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber)
                disp('Error: infinite number of iteration.')
                out = 0;
            elseif wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber<1
                disp('Error: negative number of iteration.')
                out = 0;
            end
            for ind = 1:length(wavegen.waveform.waveList)
                % test for heating
                energy = CalculateHeatEnergy(wavegen,ind);
                cooling = CalculateCoolingEnergy(wavegen,ind);
                if energy > cooling
                    coolTime = (energy - cooling)/wavegen.waveform.generator.maxAveragePower;
                    out = 0;
                    disp(['Error: maximum power density exceeded. Increase the cooling time by at least ' num2str(coolTime) '.'])
                end
                if ~isequal(length(wavegen.waveform.waveList{ind}.field),length(wavegen.waveform.waveList{ind}.time),length(wavegen.waveform.waveList{ind}.naTrigger),length(wavegen.waveform.waveList{ind}.evoTrigger))
                    out = 0;
                    disp('Error: the fields do not have the same length.')
                end
                if sum(~isfinite(wavegen.waveform.waveList{ind}.field))||sum(~isfinite(wavegen.waveform.waveList{ind}.time))
                    out = 0;
                    disp('Error: non-finite values found.')
                end
                if wavegen.waveform.waveList{ind}.field(end) ~= 0
                    out = 0;
                    disp('Error: the field does not finish at 0 T.')
                end
                dBdt = abs(diff(wavegen.waveform.waveList{ind}.field)./diff(wavegen.waveform.waveList{ind}.time));
                if sum(dBdt>wavegen.dBdtMax)
                    out = 0;
                    disp('Error: maximum slope exceeded.')
                end
                if sum(wavegen.waveform.waveList{ind}.field>wavegen.maxField)
                    out = 0;
                    disp('Error: maximum field exceeded.')
                end
                if sum(wavegen.waveform.waveList{ind}.field<wavegen.minField)
                    out = 0;
                    disp('Error: minimum field exceeded.')
                end
                if length(unique(wavegen.waveform.waveList{ind}.time))<length(wavegen.waveform.waveList{ind}.time)
                    out = 0;
                    disp('Error: some time points are repeated.')
                end
            end
        end
    end
end

