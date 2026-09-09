classdef FingerprintingGenerator < WaveformGenerator
    % this generator class can be used to generate pre-polarised and
    % non-polarised experiments with multiple evolution times, but without
    % inversion
    
    properties
        alphaArray              % list of random flip angles
        durationLimit = 1.5;    % maximum duration of an acquisition block, for heating control
        repetitionTime          % array of repetition times
        linesPerBloc = 64;      % acquisitions per bloc
    end
    
    methods
        
        function wavegen = FingerprintingGenerator(varargin)
            wavegen@WaveformGenerator;
            wavegen.name = 'Fingerprinting waveforms';
            wavegen.repetitionTime = 0.1;
        end
        
        % generate the waveforms
        function wavegen = MakeWaveform(wavegen)
            % This is the basic waveform used by default for waveform
            % generators. It is a simple polarisation - evolution
            % - detection pattern using the parameters provided. 
                                    
            % reinitialise the fields
            wavegen.waveformBlank.waveList = {};
            wavegen.waveformIndex = 0;
                
            % find the number of experiments
            pulseNumber = GetPprParameter(wavegen.waveform.pulse.pprParamList,'no_views',1);
            if length(wavegen.repetitionTime)~=pulseNumber
                wavegen.repetitionTime = wavegen.repetitionTime(1) * ones(1,pulseNumber);
            end
            
            % generate the waveforms until the heating becomes close to
            % maximum, then make a pause and start the next block
            indK = 0;
            wavegen.waveformIndex = 0;
            while indK< length(pulseNumber)
                % start a new block
                % select the next index and create the list
                wavegen.waveformIndex = wavegen.waveformIndex +1;
                wavegen = InitialiseWaveform(wavegen);
                wavegen = Delatch(wavegen);
                % polarisation period:
                if (wavegen.Tpol~=0)
                    wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bread);  % go to polarisation 
                    wavegen = Plateau(wavegen,wavegen.Tpol);      % stay at the polarisation field for the polarisation time
                end
                % go to the readout field
                wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bread);
                % fingerprinting block
                while ((wavegen.waveform.waveList{wavegen.waveformIndex}.time(end) + wavegen.repetitionTime(1)) < wavegen.durationLimit)&&(indK<pulseNumber) 
                    indK = indK + 1;
                    wavegen = PlateauWithPulse(wavegen,wavegen.Tread,wavegen.prepulseDelay,0);  % repetition time

                end
                % ends the block
                wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
                if wavegen.Tcool < minCoolTime(wavegen)
                    wavegen = Plateau(wavegen,minCoolTime(wavegen)*1.05);
                else
                    wavegen = Plateau(wavegen,wavegen.Tcool);
                end
                wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber = 1; % each bloc only plays once
            end
            
            % find out the number of experiments to program the EVO. This
            % has to be done before the undersampling.
            wavegen.waveform.experimentNumber = indK;
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
        
        function wavegen = UnderSample(wavegen)
            
            % No undersampling methods available yet for fingerprinting
            
        end
        
    end
    
end