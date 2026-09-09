classdef WaveformProfile < handle
    %WaveformProfile creates an waveformect that facilitates the creation of
    %field profiles for FFC-MRI
    %   Detailed explanation goes here
    
    properties
        waveList = {};              % list of evolution fields and times
        averageNumber = 1;          % number of repetitions for each profile
        experimentNumber = 1;       % number of experiments (number given 
                                    % to the EVO console, may not be )
                                    % the number of iterations in the case
                                    % of partial k-space acquisition)
        generator                   % fuction used to generate the waveforms
        pulse                       % contains a reference to the pulse object
        niComputerIP = '10.0.0.50'; % IP address of the NI computer (local network)
        WaveTransmitPort = 2059;    % transmission port for the waveforms.
    end
    properties (Access = private)
        figureHandle = figure('WindowStyle','Docked','Visible','off');
    end
    
    methods
        function waveform = WaveformProfile(varargin)
            % varargin = {waveformIndex timeList fieldList evoTrigger dBdt iterationNumber naVector
            % naTrigger naDelay}
            % this constructor method populates the waveform field. If the
            % object does not exist, it creates one.
            if nargin == 0
                waveform.generator = WaveformGenerator;
                waveform.generator.waveform = waveform;  % pass a reference of this object to easily access system parameters
                waveform.generator.Tevo = 0;    % by default, the waveform profile is a simple polarisation acquisition.
                waveform.generator.T1est = 0;
            else
                waveform.generator = varargin{1};
            end
        end
        
        % generate a copy of the waveform object that does not link to the
        % original (because of the handle class)
        function wavecopy = Copy(waveform)
            wavecopy = WaveformProfile;
            fieldList = fieldnames(waveform);
            for ind = 1:length(fieldList)
                wavecopy = setfield(wavecopy,fieldList{ind},getfield(waveform,fieldList{ind})); %#ok<*SFLD,*GFLD>
            end
%             wavecopy.pulse = []; % make sure the pulse object is not copied, as it is a handle
        end
        
        % transfer a copy to an existing object
        function waveform = TransferWave(waveform,wavecopy)
            fieldList = fieldnames(waveform);
            for ind = 1:length(fieldList)
                waveform = setfield(waveform,fieldList{ind},getfield(wavecopy,fieldList{ind})); %#ok<*SFLD,*GFLD>
            end
        end
        
        % test for safety and runs the generator to make the waveforms
        function out = MakeWaveform(waveform)
            waveform.generator = Update(waveform.generator);
            if ~TestSequence(waveform.generator)
                disp('Warning: sequence unsafe to be played.')
                out = 0;
            else
                out = 1;
            end
        end
                
        function time = FindOverallDuration(waveform)
            time = 0;
            for ind = 1:length(waveform.waveList)
                time = time + waveform.waveList{ind}.time(end);
            end
        end
                
        function plot(waveform,varargin)
            if ~ishandle(waveform.figureHandle)
                waveform.figureHandle = figure('WindowStyle','Docked');
            end
            set(waveform.figureHandle,'Visible','on')
            figure(waveform.figureHandle)
            if isempty(varargin)
                for ind = 1:length(waveform.waveList)
                    plot(waveform.waveList{ind}.time,waveform.waveList{ind}.field,waveform.waveList{ind}.time,waveform.waveList{ind}.naTrigger/10.01,waveform.waveList{ind}.time,waveform.waveList{ind}.evoTrigger/10)
                    hold on
                    if isfield(waveform.waveList{ind},'LETIprotect')
                        plot(waveform.waveList{ind}.time,waveform.waveList{ind}.LETIprotect/10.2,'k')
                    end
                end
                hold off
            else
                ind = varargin{1};
                plot(waveform.waveList{ind}.time,waveform.waveList{ind}.field,waveform.waveList{ind}.time,waveform.waveList{ind}.naTrigger/5.01,waveform.waveList{ind}.time,waveform.waveList{ind}.evoTrigger/5)
            end
        end
        
        function out = UploadWaveform(waveform)
            % out = UploadWaveform(waveform)
            % Transmit the waveform to the NI computer
            if NI_client(waveform) %.niComputerIP, waveform.WaveTransmitPort, waveform.averageNumber, waveform.iterationNumber, waveform.naVector,...
%                              waveform.time, waveform.field, waveform.naTrigger, waveform.evoTrigger);
                out = 1;
                disp('Waveform transmitted successfully.')
            else
                out = 0;
                disp('Error during the transmission of the waveform.')
            end
        end
        
    end    
    
end

