classdef PPLookLockerGenerator < WaveformGenerator
    % this generator class can be used to generate pre-polarised and
    % non-polarised experiments with multiple evolution times, but without
    % inversion
    
    properties
        no_evolution_times %number of times to loop over
    end
    
    methods
        
        function wavegen = PPLookLockerGenerator(varargin)
            wavegen@WaveformGenerator;
            wavegen.name = 'Prepolarised Look Locker';
            wavegen.Tpol = 0.4; % polarised sequence by default
            wavegen.maxAveragePower = 40e3;
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
            wavegen.Tevo = wavegen.Tread+wavegen.readDelay;
%             wavegen.Tevo = wavegen.Tevo*ones(1,GetPprParameter(wavegen.waveform.pulse.pprParamList,'no_echoes',2)-1);
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
            wavegen.waveformBlank.waveList = {};
            wavegen.waveformIndex = 0;
            wavegen.waveform.waveList = [];
                
            % estimation of the relaxation time at readout field
           
            
            if ((wavegen.Tread+wavegen.readDelay)*(wavegen.no_evolution_times)) > 1.5
                error('Sequence too harsh, the magnet will melt.')
            end
            if wavegen.Tread > wavegen.Tread
                disp(['Additional time required:' num2str(wavegen.Tread - wavegen.Tread) 's'])
                error('Acquisition time too short, some slices may not be recorded correctly')
            end
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
                    wavegen.dBdtPerf = 10;
                    wavegen.readDelay = 0.020;
                    % Start the waveform shape here:
                    wavegen = InitialiseWaveform(wavegen);
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Bevo = wavegen.Bevo(indField);
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Tevo = wavegen.Tevo(indField,indTevo);
                    
                    % make the waveform segment by segment
                    wavegen = Delatch(wavegen);   % adds the delatch signal for the IECO amps
                    % polarisation period:
                    if (wavegen.Tpol~=0)
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bpol);  % go to polarisation
                        wavegen = Plateau(wavegen,wavegen.Tpol);      % stay at the polarisation field for the polarisation time
                    end
                    for n=1:wavegen.no_evolution_times-1
                        % now going to evolution field:
                        if wavegen.Tevo(indField,indTevo) ~= 0
                            wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bevo(indField));
                            wavegen = Plateau(wavegen,wavegen.Tevo(indField,indTevo));
                        end
                        % and finally go to readout
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bread);
                        wavegen = PlateauWithPulse(wavegen,wavegen.Tread,wavegen.readDelay,0);
                        
                        % add the navigator for the last scan
                        if n == wavegen.no_evolution_times-1
                            wavegen = PlateauWithPulse(wavegen,wavegen.Tread,0.01,0);
                        end
                            
                    end
                    % ends with the cooling time
                    wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
                    if wavegen.Tcool < minCoolTime(wavegen)
                        wavegen = Plateau(wavegen,minCoolTime(wavegen)*1.1);
                    else
                        wavegen = Plateau(wavegen,wavegen.Tcool);
                    end
                    wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber = GetPprParameter(wavegen.waveform.pulse.pprParamList,'no_views',1);
                end
                
            end
            % find out the number of experiments to program the EVO. This
            % has to be done before the undersampling.
            wavegen.waveform.experimentNumber = wavegen.waveformIndex;
            wavegen.waveform.averageNumber = wavegen.averageNumber;
            
        end
        
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
    end
    
end