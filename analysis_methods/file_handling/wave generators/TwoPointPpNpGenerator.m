classdef TwoPointPpNpGenerator < WaveformGenerator
    % this generator class can be used to generate two-points experiment
    % profiles
    
    properties
        inversionList = [];     % list of booleans describing the presence of inversion for each Bevo
        Tinv = 0.001;           % duration of the inversion pulse
        prepolarised = [];      % list of prepolarised experiments
    end
    
    methods
        
        function wavegen = TwoPointPpNpGenerator(varargin)
            wavegen@WaveformGenerator;
            wavegen.name = '2-pts PP/NP';
        end
        
        % For the two-pt method we only need one time point per evolution
        % field
        function wavegen = UpdateTime(wavegen)
            wavegen.Tevo = round(wavegen.T1est*1e6)/1e6; % limits the precision to 1us
        end
        
         % generate the waveforms for 2-pts experiments:
         % - first make a non-inverted measurement at polarisation field
         % (and detection field if they differ)
         % - then do inverted measurement at these fields, plus all the
         % others
        function wavegen = MakeWaveform(wavegen)
                                    
            % reinitialise the fields
            wavegen.waveformBlank.waveList = {};
            wavegen.waveformIndex = 0;
            
            % list the pre-polarised experiments
            wavegen.prepolarised = wavegen.Bevo < wavegen.Bpol/2;
            
            % add the pre-polarised measurements at detection and
            % polarisation fields
            if wavegen.Bpol == wavegen.Bread
                wavegen.prepolarised = [1;wavegen.prepolarised(:)];
                wavegen.Bevo = [wavegen.Bpol;wavegen.Bevo(:)];
            else
                wavegen.prepolarised = [1;1;wavegen.prepolarised(:)];
                wavegen.Bevo = [wavegen.Bpol;wavegen.Bread;wavegen.Bevo(:)];
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
                    
                    % Start the waveform shape here:
                    wavegen = InitialiseWaveform(wavegen);                    
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Bevo = wavegen.Bevo(indField);
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Tevo = wavegen.Tevo(indField,indTevo);
                    
                    % make the waveform segment by segment
                    wavegen = Delatch(wavegen);   % adds the delatch signal for the IECO amps
                    % polarisation period:
                    if (wavegen.Tpol~=0)&&(wavegen.prepolarised(indField))
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bpol);  % go to polarisation 
                        wavegen = Plateau(wavegen,wavegen.Tpol);      % stay at the polarisation field for the polarisation time
                    end
                    % inversion (if required)
                    if wavegen.inversionList(indField)
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bevo(indField));
                        wavegen = PlateauWithPulse(wavegen,wavegen.Tinv,wavegen.readDelay,0);
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

                end
            end
            % find out the number of experiments to program the EVO. This
            % has to be done before the undersampling.
            wavegen.waveform.experimentNumber = wavegen.waveformIndex;
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
            SetPprData(wavegen.waveform.pulse,'inv_obs_mod_level', 1);
        end
    end
    
end