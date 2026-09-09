classdef ConstantField < WaveformGenerator
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
        BDClim = 0.12;
    end
    
    methods
        function wavegen = ConstantField(varargin)
            wavegen@WaveformGenerator;
            wavegen.name = 'DC field';
            wavegen.Bread = 0.1;
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
            
            % check that the field is below 120 mT
            wavegen.Bread = GetFrequencyInHz(wavegen.waveform.pulse.pprParamList)/gammaH; % field in T
            if abs(wavegen.Bread) > 0.12
                error('Field too high for DC mode')
            end
            
            % switch on the field
            wavegen.waveformIndex = 1;
            wavegen = InitialiseWaveform(wavegen); 
            wavegen = Delatch(wavegen);
            wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,wavegen.Bread);
%             
%             % make all the acquisitions
%             wavegen.waveformIndex = 2;
%             wavegen = InitialiseWaveform(wavegen);
%             wavegen.waveform.waveList{2}.time = [0 2];
%             wavegen.waveform.waveList{2}.field = [wavegen.Bread wavegen.Bread];
%             wavegen.waveform.waveList{2}.naTrigger = [0 0];
%             wavegen.waveform.waveList{2}.evoTrigger = [0 0];
%             wavegen.waveform.waveList{2}.iterationNumber = 1;
%             wavegen.waveform.waveList{2}.naVector = [0 0 0];
%             wavegen.waveform.waveList{2}.Bevo = wavegen.Bread;
%             wavegen.waveform.waveList{2}.Tevo = 2;
%             wavegen.waveform.waveList{2}.averageNumber = 1;
%             wavegen.waveform.waveList{2}.indexBevo = 1;
%             wavegen.waveform.waveList{2}.indexTevo = 1;
            
            % switch off the field
            wavegen.waveformIndex = 2;
            wavegen = InitialiseWaveform(wavegen);
            
            wavegen = InitialiseWaveform(wavegen);
            wavegen.waveform.waveList{wavegen.waveformIndex}.time = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.field = wavegen.Bread;
            wavegen.waveform.waveList{wavegen.waveformIndex}.naTrigger = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.evoTrigger = 0;
            wavegen.waveform.waveList{wavegen.waveformIndex}.iterationNumber = 1;
            wavegen.waveform.waveList{wavegen.waveformIndex}.naVector = [0 0 0];
            wavegen.waveform.waveList{wavegen.waveformIndex}.Bevo = wavegen.Bread;
            wavegen.waveform.waveList{wavegen.waveformIndex}.Tevo = 2;
            wavegen.waveform.waveList{wavegen.waveformIndex}.averageNumber = 1;
            wavegen.waveform.waveList{wavegen.waveformIndex}.indexBevo = 1;
            wavegen.waveform.waveList{wavegen.waveformIndex}.indexTevo = 1;
            wavegen = RampConstantSlope(wavegen,wavegen.dBdtComfort,0);
              
        end        
        
        % undersample the k-space. kSpaceLine is an array of booleans that
        % is set to 'true' if the line is measured, 'false' if it is
        % skipped
        function wavegen = UnderSample(wavegen)
            
        end
        
        % this function holds custom scripts that run after the generator
        % has finished creating the waveforms, it can be used to update
        % some fields in the PPL file.
        function wavegen = UpdatePulse(wavegen)
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Safety testing functions

        
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
                if abs(wavegen.waveform.waveList{ind}.field(end)) > wavegen.BDClim
                    out = 0;
                    disp('Error: the field does not finish below 0.12 T.')
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

