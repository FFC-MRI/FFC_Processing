classdef InversionRecoveryGeneratorSoftRamp < InversionRecoveryGenerator
    % this generator class can be used to generate inversion recovery
    % experiments at different evolution times
    
    properties 
        dBdtSoft = 1;
        softRampDuration = 1e-3;
    end
    
    methods
        
        function wavegen = InversionRecoveryGeneratorSoftRamp(varargin)
            wavegen@InversionRecoveryGenerator;
            wavegen.name = 'Inversion Recovery with soft ramps';
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
                    
                    % make the waveform segment by segment
                    wavegen = Delatch(wavegen);   % adds the delatch signal for the IECO amps
                    % polarisation period:
                    if wavegen.Tpol~=0
%                         wavegen = BrokenRamp(wavegen,wavegen.dBdtComfort,wavegen.dBdtSoft,wavegen.softRampDuration,wavegen.Bpol);
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bpol);  % go to polarisation 
                        wavegen = Plateau(wavegen,wavegen.Tpol);      % stay at the polarisation field for the polarisation time
                    end
                    % inversion
                    if (wavegen.Tinv~=0)&&(wavegen.inversionFlag(indField))
%                         wavegen = BrokenRamp(wavegen,wavegen.dBdtPerf,wavegen.dBdtSoft,wavegen.softRampDuration,wavegen.Bread);
                        wavegen = RoundedRamp(wavegen,wavegen.dBdtPerf,wavegen.softRampDuration,wavegen.Bread);
                        wavegen = PlateauWithPulse(wavegen,wavegen.Tinv,wavegen.readDelay,0.001);
                    end
                    % now going to evolution field:
                    if wavegen.Tevo(indField,indTevo) ~= 0
%                         wavegen = BrokenRamp(wavegen,wavegen.dBdtPerf,wavegen.dBdtSoft,wavegen.softRampDuration,wavegen.Bevo(indField));
                        wavegen = RoundedRamp(wavegen,wavegen.dBdtPerf,wavegen.softRampDuration,wavegen.Bevo(indField));
                        wavegen = Plateau(wavegen,wavegen.Tevo(indField,indTevo));
                    end
                    % and finally go to readout                    
%                     wavegen = BrokenRamp(wavegen,wavegen.dBdtPerf,wavegen.dBdtSoft,wavegen.softRampDuration,wavegen.Bread);
                    wavegen = RoundedRamp(wavegen,wavegen.dBdtPerf,wavegen.softRampDuration,wavegen.Bread);
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
        
        
        % makes a ramp with a slower ending segment to reduce the voltage
        % and ringing time
        function wavegen = BrokenRamp(wavegen,slope1,slope2,time2,fieldEnd)
            % find the field that has to be reached between the two
            % segments
            Bini = wavegen.waveform.waveList{wavegen.waveformIndex}.field(end);
            dB = fieldEnd - Bini;
            time1 = abs(abs(dB) - slope2*time2)/slope1;
            fieldInter = sign(fieldEnd-Bini)*time1*slope1 + Bini;
            wavegen = RampConstantSlope(wavegen,slope1,fieldInter);
            wavegen = RampConstantSlope(wavegen,slope2,fieldEnd);
        end
        
        % Same as above but the limit between the two ramps is provided by
        % the value of the magnetic field instead of the time
        function wavegen = BrokenRampField(wavegen,slope1,slope2,field2,fieldEnd)
            % find the field that has to be reached between the two
            % segments
            wavegen = RampConstantSlope(wavegen,slope1,field2);
            wavegen = RampConstantSlope(wavegen,slope2,fieldEnd);
        end
        
        % same as the two above, but using rounded edges using a circle,
        % with continuous derivative
        function wavegen = RoundedRamp(wavegen,slope,roundTime,fieldEnd)
            dt = 0.1e-3; % resolution of the rounded element
            time = 0:dt:roundTime;
            Bini = wavegen.waveform.waveList{wavegen.waveformIndex}.field(end);
            if Bini == fieldEnd
                return
            end
            [B,Bdt] = FindSoftCornerParameter(slope,roundTime,time,0);
            if Bini < fieldEnd % positive slope
                
            else % negative slope
                B = -B;
                Bdt = -Bdt;
            end
            % make sure that there is enough time to make a rounded edge
            if abs(Bdt)>abs(Bini - fieldEnd)
                wavegen = RampConstantSlope(wavegen,slope,fieldEnd);
                return
            end
            % add the rounded edge
            for pt = 1:length(B)
                wavegen.waveform.waveList{wavegen.waveformIndex}.time(end+1) = wavegen.waveform.waveList{wavegen.waveformIndex}.time(end) + dt;
                wavegen.waveform.waveList{wavegen.waveformIndex}.field(end+1) = B(pt) + Bini;
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(end+1) = 0;
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(end+1) = 0;
            end
            % continue the ramp normally
            wavegen = RampConstantSlope(wavegen,slope,fieldEnd - Bdt);
            % add the rounded edge
            for pt = 2:length(B)
                wavegen.waveform.waveList{wavegen.waveformIndex}.time(end+1) = wavegen.waveform.waveList{wavegen.waveformIndex}.time(end) + dt;
                wavegen.waveform.waveList{wavegen.waveformIndex}.field(end+1) = fieldEnd - B(end-pt+1);
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(end+1) = 0;
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(end+1) = 0;
            end
        end
        
        
    end
    
end