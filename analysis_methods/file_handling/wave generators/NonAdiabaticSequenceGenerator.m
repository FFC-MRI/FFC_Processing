classdef NonAdiabaticSequenceGenerator < WaveformGenerator
    % this generator class can be used to generate two-points experiment
    % profiles
    
    
    properties      
    end
    
    methods
                
        function wavegen = NonAdiabaticSequenceGenerator(varargin)
            wavegen@WaveformGenerator;
            wavegen.name = 'Inversion Recovery with non-adiabatic transition control';
        end
        
        function wavegen = SetNonAdiabaticVector(wavegen,vector)
            % test the inputs
            
            % assign the value
            wavegen.waveform.waveList{wavegen.waveformIndex}.naVector = vector;
            
        end
        
        function wavegen = PlateauWithNonAdiabaticField(wavegen,time)
            B = wavegen.waveform.waveList{wavegen.waveformIndex}.field(end);
            % add the ramp to the current waveform object
            timeStamp =                                                   [                                                              time+wavegen.minDeltaT   wavegen.minDeltaT];
            wavegen.waveform.waveList{wavegen.waveformIndex}.field =      [wavegen.waveform.waveList{wavegen.waveformIndex}.field        B      B];
            wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger =  [wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger    1      1];
            wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger = [wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  0       0];
             wavegen.waveform.waveList{wavegen.waveformIndex}.time =      [wavegen.waveform.waveList{wavegen.waveformIndex}.time        wavegen.waveform.waveList{wavegen.waveformIndex}.time(end)+cumsum(timeStamp)];
        end
        
        function wavegen = RampNonAdiabatic(wavegen,slope,delay,fieldEnd)
            % waveform = RampConstantTime(waveform,rampDuration,delay,fieldEnd)
            % This method appends a slope with a non-adiabatic transition
            % when reaching the lowest field. A delay can be added to
            % maintain the transitory field a little longer.
            % this must be programmed to follow adiabatic field profiles!!!
            
            % test the inputs
            if ~isobject(wavegen)
                error('First input must be a WaveformGenerator object')
            end
            if (slope > wavegen.dBdtMax)||(slope < wavegen.dBdtMin)
                error('The field time-derivate is out of boundaries.')
            end
            if (fieldEnd > wavegen.maxField)||(fieldEnd < wavegen.minField)
                error('The field is out of boundaries.')
            end
            
            rampDuration = abs(wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)-fieldEnd)/slope;
            % find if the transition is at the start or at the end
            if fieldEnd > wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)    % transition is at the start
                % add the ramp to the current waveform object
                timeStamp = [                     wavegen.minDeltaT  delay+wavegen.minDeltaT           rampDuration  wavegen.minDeltaT];
                wavegen.waveform.waveList{wavegen.waveformIndex}.field = [wavegen.waveform.waveList{wavegen.waveformIndex}.field            wavegen.waveform.waveList{wavegen.waveformIndex}.field(end) wavegen.waveform.waveList{wavegen.waveformIndex}.field(end)  fieldEnd      fieldEnd];
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger = [wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger    1              1               1             0];
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger = [wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  0              0               0             0];
            else                            % transition is at the end
                timeStamp = [                     wavegen.minDeltaT  rampDuration   delay+wavegen.minDeltaT     wavegen.minDeltaT];
                wavegen.waveform.waveList{wavegen.waveformIndex}.field = [wavegen.waveform.waveList{wavegen.waveformIndex}.field            wavegen.waveform.waveList{wavegen.waveformIndex}.field(end) fieldEnd       fieldEnd  fieldEnd];
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger = [wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger    1              1              1         1];
                wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger = [wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  0              0              0         0];
            end
            wavegen.waveform.waveList{wavegen.waveformIndex}.time = [wavegen.waveform.waveList{wavegen.waveformIndex}.time wavegen.waveform.waveList{wavegen.waveformIndex}.time(end)+cumsum(timeStamp)];
        end
        
        
    end
    
end
