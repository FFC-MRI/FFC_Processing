classdef DoubleResonanceGenerator < WaveformGenerator
    % this generator class can be used to generate two-points experiment
    % profiles
    
    
    properties      
    end
    
    methods
        
           
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
                        wavegen = PlateauWithPulse(wavegen,wavegen.Tevo(indField,indTevo),0,0);
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
        
        
        function waveform = AcquistionWithFrequencyShift(waveform,pulseDuration,acqDuration,acquitionDelay,finishDelay,Btrans,Bread)
            % Generates a plateau with two segments, for acquisition at
            % frequencies different from the excitation. This allows
            % avoiding ringing.
            % test the inputs
            
            waveform = RampConstantSlope(waveform,waveform.dBdtPerf,Btrans);
            indIni = length(waveform.waveList{waveform.waveformIndex}.evoTrigger);
            waveform = PlateauWithPulse(waveform,pulseDuration,acquitionDelay,0);
            waveform = RampConstantSlope(waveform,waveform.dBdtMax*0.95,Bread);
            waveform = PlateauWithPulse(waveform,acqDuration,0,finishDelay);
            waveform.waveList{waveform.waveformIndex}.evoTrigger((indIni+3):(end-2)) = 1;
            
        end
    end
end