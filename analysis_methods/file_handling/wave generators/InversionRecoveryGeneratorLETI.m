classdef InversionRecoveryGeneratorLETI < InversionRecoveryGenerator & NonAdiabaticSequenceGenerator
    % this generator class can be used to generate inversion recovery
    % experiments at different evolution times
    
    properties 
        dBdtSoft = 3; % 2.5 A/ms
        Blimit = 1.224e-3; %1e-4; % lower limit of the field below which the LETI PSU is enabled
        LETISettleTime = 5e-3; % settling time of the LETI PSU after enabling
        LETIAdditionalDelay = 1e-3; % additional delay to allow the IECO voltage to reduce when starting a broken ramp
        LETIProtectDelay = 20e-3; % delay to escape protection mode
        LETIScalingFactor = 10/1e-4;      % scaling factor for the LETI PSU (in V/T)
        LETIprotect = [];
    end
    
    methods
        
        function wavegen = InversionRecoveryGeneratorLETI(varargin)
            wavegen@InversionRecoveryGenerator;
            wavegen.name = 'Inversion Recovery for LETI PSU';
        end
                
        % generate the waveforms
        function wavegen = MakeWaveform(wavegen)
            % This is the basic waveform used by default for waveform
            % generators. It is a simple polarisation - evolution
            % - detection pattern using the parameters provided. 
                                    
            % reinitialise the fields
            wavegen.waveformBlank.waveList = {};
            wavegen.waveformIndex = 0;
            wavegen.Blimit = 1.3e-3;
                            
            % generate the waveforms
            for indField = 1:length(wavegen.Bevo)
                for indTevo = 1:size(wavegen.Tevo,2)
                    % select the next index and create the list
                    wavegen.waveformIndex = wavegen.waveformIndex +1;
                    ind = wavegen.waveformIndex;
                    
                    % Start the waveform shape here:
                    wavegen = InitialiseWaveform(wavegen);                    
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Bevo = wavegen.Bevo(indField);
                    wavegen.waveform.waveList{wavegen.waveformIndex}.Tevo = wavegen.Tevo(indField,indTevo);
                    wavegen.waveform.waveList{wavegen.waveformIndex}.averageNumber = wavegen.averageNumber;
                    % store the index of the evolution field and time into
                    % an array:
                    wavegen.waveform.waveList{ind}.indexBevo = indField;
                    wavegen.waveform.waveList{ind}.indexTevo = indTevo;
                    wavegen.indexList(ind,:) = [indField indTevo];
                    
                    % finds when the LETI PSU is activated
                    if wavegen.Bevo(indField) > wavegen.Blimit
                        field = wavegen.Bevo(indField);
                         wavegen.waveform.waveList{ind}.naVector = [0 0 0];
                    else
                        field = 0;
                        wavegen.waveform.waveList{ind}.naVector = wavegen.LETIScalingFactor*[0 0 wavegen.Bevo(indField)];
                    end
                    
                    % make the waveform segment by segment
                    wavegen = Delatch(wavegen);   % adds the delatch signal for the IECO amps
                    % polarisation period:
                    if wavegen.Tpol~=0
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bpol);  % go to polarisation 
                        wavegen = Plateau(wavegen,wavegen.Tpol);      % stay at the polarisation field for the polarisation time
                    end
                    % inversion
                    if (wavegen.Tinv~=0)&&(wavegen.inversionFlag(indField))
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bread);
                        wavegen = PlateauWithPulse(wavegen,wavegen.Tinv,wavegen.readDelay,wavegen.inversionTime);
                    end
                    % now going to evolution field:
                    if wavegen.Bevo(indField)>wavegen.Blimit  % case when the IECO amplifiers can make the evolution field
                        if wavegen.Tevo(indField,indTevo) ~= 0
                            wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,field);
                            wavegen = Plateau(wavegen,wavegen.Tevo(indField,indTevo));
                        end
                        wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bread);
                    else   % case when the evolution field is maintained by the LETI PSU
                        if wavegen.Tevo(indField,indTevo) ~= 0
                            wavegen = BrokenRampLeti(wavegen,wavegen.dBdtPerf,wavegen.dBdtSoft,wavegen.LETISettleTime,0);
                            wavegen = PlateauWithNonAdiabaticField(wavegen,wavegen.Tevo(indField,indTevo));
                            wavegen.waveform.waveList{wavegen.waveformIndex}.naVector(3) = wavegen.Bevo(indField);
                        end
                        wavegen = BrokenRampLeti(wavegen,wavegen.dBdtPerf,wavegen.dBdtSoft,wavegen.LETISettleTime,wavegen.Bread);
%                         wavegen = RampConstantSlope(wavegen,wavegen.dBdtPerf,wavegen.Bread);
                    end
                    % and finally go to readout 
                    wavegen = PlateauWithPulse(wavegen,wavegen.Tread,wavegen.readDelay,0);
                    % remove duplicate time points
                    wavegen = removeduplicate(wavegen);
                    % ends with the cooling time
                    wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
                    if wavegen.Tcool < minCoolTime(wavegen)
                        wavegen = Plateau(wavegen,minCoolTime(wavegen)*1.05);
                    else
                        wavegen = Plateau(wavegen,wavegen.Tcool);
                    end
                    if wavegen.Bevo(indField)<wavegen.Blimit
                        wavegen = GenerateProtectLine(wavegen);
                        wavegen = removeduplicate(wavegen);
                    end
                    
                end
            end
            % find out the number of experiments to program the EVO. This
            % has to be done before the undersampling.
            wavegen.waveform.experimentNumber = wavegen.waveformIndex;
            wavegen.waveform.averageNumber = wavegen.averageNumber;
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
        
        % makes a ramp with a slower ending segment to reduce the voltage
        % and ringing time
        function wavegen = BrokenRampLeti(wavegen,slope1,slope2,time2,fieldEnd)
            % find the field that has to be reached between the two
            % segments            
            Bini = wavegen.waveform.waveList{wavegen.waveformIndex}.field(end);
            dB = fieldEnd - Bini;
            if Bini>fieldEnd
                time1 = abs(abs(dB) - slope2*time2)/slope1;
                fieldInter = sign(fieldEnd-Bini)*time1*slope1 + Bini;
                if fieldInter>fieldEnd
                    wavegen = RampConstantSlope(wavegen,slope1,fieldInter);
                    wavegen = RampNonAdiabatic(wavegen,slope2,0,fieldEnd);
                else
                    wavegen = RampNonAdiabatic(wavegen,slope2,0,fieldEnd);
                end
            else
                fieldInter = sign(fieldEnd-Bini)*time2*slope2 + Bini;
                if fieldInter<fieldEnd
                    wavegen = RampNonAdiabatic(wavegen,slope2,0,Bini+1e-9);
%                     wavegen = RampConstantSlope(wavegen,slope2,fieldInter);
                    wavegen = RampConstantSlope(wavegen,slope1,fieldEnd);
                else
                    wavegen = RampConstantSlope(wavegen,slope1,fieldEnd);
                end
            end
        end
        
        % Same as above but the limit between the two ramps is provided by
        % the value of the magnetic field instead of the time
        function wavegen = BrokenRampField(wavegen,slope1,slope2,field2,fieldEnd)
            % find the field that has to be reached between the two
            % segments
            wavegen = RampConstantSlope(wavegen,slope1,field2);
            wavegen = RampConstantSlope(wavegen,slope2,fieldEnd);
        end
        
        function wavegen = GenerateProtectLine(wavegen)
          
            % add the delay for LETI PSU
            ind = find(wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger);
            timestart = wavegen.waveform.waveList{wavegen.waveformIndex}.time(ind(1));
            timedelayed = timestart + wavegen.LETIAdditionalDelay + [0 1e-6];
            if wavegen.LETIAdditionalDelay>0
                val = interp1(wavegen.waveform.waveList{wavegen.waveformIndex}.time,wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger,timedelayed(end));
            else
                val = 1;
            end
            wavegen.waveform.waveList{wavegen.waveformIndex}.field       = [wavegen.waveform.waveList{wavegen.waveformIndex}.field,      interp1(wavegen.waveform.waveList{wavegen.waveformIndex}.time,wavegen.waveform.waveList{wavegen.waveformIndex}.field,timedelayed)];
            wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  = [wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger, interp1(wavegen.waveform.waveList{wavegen.waveformIndex}.time,wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger,timedelayed)];
%             wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect = [wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect,interp1(wavegen.waveform.waveList{wavegen.waveformIndex}.time,wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect,timedelayed)];
            wavegen.waveform.waveList{wavegen.waveformIndex}.time        = [wavegen.waveform.waveList{wavegen.waveformIndex}.time,       timedelayed];
            wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger   = [wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger,  [0 val]];
            [wavegen,indexsorted] = sortarray(wavegen);
            if ind(1)<indexsorted(end-1)
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(ind(1):indexsorted(end-1)) = 0;
            else
                wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(indexsorted(end-1)-1:ind(1)+1) = 1;
            end
            
            
            % generate the LO/HI VHT trigger
            ind = find(wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger);
            timestartnew = wavegen.waveform.waveList{wavegen.waveformIndex}.time(ind(1));
            
            timestart = max(timestart,timestartnew); % the LO/HI VHT signal is only needed 20 ms before the end of the broken ramp anyway
            
%             timeenable = timestart - wavegen.LETIProtectDelay;
            timeenable = timestart - wavegen.LETIProtectDelay + [0,1e-6, 10e-3, 10.001e-3];
            i = wavegen.waveformIndex;
            
            wavegen.waveform.waveList{i}.field       = [wavegen.waveform.waveList{i}.field,      interp1(wavegen.waveform.waveList{i}.time,wavegen.waveform.waveList{i}.field,timeenable)];
            wavegen.waveform.waveList{i}.naTrigger   = [wavegen.waveform.waveList{i}.naTrigger,  interp1(wavegen.waveform.waveList{i}.time,wavegen.waveform.waveList{i}.naTrigger,timeenable)];
            wavegen.waveform.waveList{i}.evoTrigger  = [wavegen.waveform.waveList{i}.evoTrigger, interp1(wavegen.waveform.waveList{i}.time,wavegen.waveform.waveList{i}.evoTrigger,timeenable)];
            wavegen.waveform.waveList{i}.time        = [wavegen.waveform.waveList{i}.time,       timeenable];
            wavegen.waveform.waveList{i}.LETIprotect = zeros(1,numel(wavegen.waveform.waveList{wavegen.waveformIndex}.time));
            wavegen.waveform.waveList{i}.LETIprotect(end-2:end-1) = 1;
            [wavegen,indexsorted] = sortarray(wavegen);
            wavegen.waveform.waveList{i}.LETIprotect(indexsorted(end-2):indexsorted(end-1)) = 1;
            
            
        end
        
%         function wavegen = insertIntoWaveform(wavegen,timeenable,position)
%             timestart = wavegen.waveform.waveList{wavegen.waveformIndex}.time(position);
%             wavegen.waveform.waveList{wavegen.waveformIndex}.time        = insertIntoArray(wavegen.waveform.waveList{wavegen.waveformIndex}.time, timeenable + [0,1e-6],position);
%             wavegen.waveform.waveList{wavegen.waveformIndex}.field       = insertIntoArray(wavegen.waveform.waveList{wavegen.waveformIndex}.field,wavegen.waveform.waveList{wavegen.waveformIndex}.field(position)*[1 1],position);
%             wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger   = insertIntoArray(wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger,wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(position)*[1 1],position);
%             wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  = insertIntoArray(wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger,wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(position)*[1 1],position);
% %             if status==1
% %                 wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect = insertIntoArray(wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect,[0 1],position);
% %             else
% %                 wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect = insertIntoArray(wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect,[1 0],position);
% %             end
%         end
        
        function [wavegen,indexsorted] = sortarray(wavegen)
            [wavegen.waveform.waveList{wavegen.waveformIndex}.time,ind]  = sort(wavegen.waveform.waveList{wavegen.waveformIndex}.time);
            wavegen.waveform.waveList{wavegen.waveformIndex}.field       = wavegen.waveform.waveList{wavegen.waveformIndex}.field(ind);
            wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger   = wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(ind);
            wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  = wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(ind);
            if isfield(wavegen.waveform.waveList{wavegen.waveformIndex},'LETIprotect')
                if numel(wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect) == numel(ind)
                    wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect = wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect(ind);
                end
            end
            indexsorted = 1:numel(ind);
            indexsorted(ind) = indexsorted;
%             indstart = find(ind == max(ind)-2);
%             indstop = find(ind == max(ind)-1);
        end
        
        function  wavegen = removeduplicate(wavegen)
            [wavegen.waveform.waveList{wavegen.waveformIndex}.time,ind] = unique(wavegen.waveform.waveList{wavegen.waveformIndex}.time);
            wavegen.waveform.waveList{wavegen.waveformIndex}.field       = wavegen.waveform.waveList{wavegen.waveformIndex}.field(ind);
            wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger   = wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger(ind);
            wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger  = wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger(ind);
            if isfield(wavegen.waveform.waveList{wavegen.waveformIndex},'LETIprotect')
                wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect = wavegen.waveform.waveList{wavegen.waveformIndex}.LETIprotect(ind);
            end
        end
    end
end