classdef PpNpIrGenerator < WaveformGenerator
    % this generator class can be used to generate inversion recovery
    % experiments at different evolution times
    
    properties 
    end
    
    methods
        
        function wavegen = PpNpIrGenerator(varargin)
            wavegen@WaveformGenerator;
            wavegen.name = 'NP/PP with fake inversion pulse';
        end
        
        function wavegen = UpdateField(wavegen)
            % make the list of evolution fields
            wavegen.Bevo = logspace(log10(wavegen.BevoStart),log10(wavegen.BevoEnd),wavegen.NBevo)'; % the first dimension is set to the number of fields
            wavegen = ProbeQP(wavegen); % adds the sampling for the QP
            wavegen.inversionFlag = true(size(wavegen.Bevo));
        end
        
        % generate the waveforms
        function wavegen = MakeWaveform(wavegen)
            % This is the basic waveform used by default for waveform
            % generators. It is a simple polarisation - evolution
            % - detection pattern using the parameters provided. 
                                    
            % reinitialise the fields
            wavegen.waveformBlank.waveList = {};
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
                    wavegen.waveform.waveList{wavegen.waveformIndex}.averageNumber = wavegen.averageNumber;
                    
                    % make the waveform segment by segment
                    wavegen = Delatch(wavegen);   % adds the delatch signal for the IECO amps
                    
                    % fake inversion
                    if (wavegen.Tinv~=0)&&(wavegen.inversionFlag(indField))
                        wavegen = PlateauWithPulse(wavegen,wavegen.Tinv,wavegen.readDelay,0.001);
                    end
                    
                    % polarisation period:
                    if (wavegen.Tpol~=0)&&(wavegen.Bevo(indField)<wavegen.Bpol/2)
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
                    % ends with the cooling time0
                    wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
                    if wavegen.Tcool < minCoolTime(wavegen)
                        wavegen = Plateau(wavegen,minCoolTime(wavegen)*1.05);
                    else
                        cooltime = max(0,wavegen.Tcool-wavegen.cDaqDelay);
                        wavegen = Plateau(wavegen,wavegen.Tcool);
                    end

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
%             SetPprData(wavegen.waveform.pulse,'no_averages',wavegen.waveform.averageNumber);
            
            % copy the waveform parameters into the PPR file
            if ~isempty(wavegen.Tevo)
                times = zeros(1,50);
                times(1:length(wavegen.Tevo(:))) = wavegen.Tevo(:)';
                SetPprData(wavegen.waveform.pulse,'t_evol',  round(times(1:50)*1e3));
            end
            if ~isempty(wavegen.Bevo)
                fields = zeros(1,50);
                fields(1:length(wavegen.Bevo)) = wavegen.Bevo(:)';
                SetPprData(wavegen.waveform.pulse,'b_evol',round(fields(1:50)*1e3));
            end
            if ~isempty(wavegen.Bpol)
                val = wavegen.Bpol;
                SetPprData(wavegen.waveform.pulse,'b_pol',round(val*1e3));
            end
            if ~isempty(wavegen.Tpol)
                val = wavegen.Tpol;
                SetPprData(wavegen.waveform.pulse,'t_pol',round(val*1000));
            end
            if ~isempty(wavegen.Tcool)
                val = wavegen.Tcool;
                SetPprData(wavegen.waveform.pulse,'down_time',round(val*1000));
            end
            % by default, no inversion
            SetPprData(wavegen.waveform.pulse,'inv_obs_mod_level', 1);
        end
        
    end
    
end