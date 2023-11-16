classdef (Abstract) R1Template

%Give the Name that will be compared to in the GUI and Checker
    properties (Abstract)
    processName;
    end


%Example Code
    methods (Abstract)
        [R1out,others,tres,sse,err,fitres] = FitRelaxation(imageSmoothed,t,B0,B0_pol,rois)
            %Insert Desired R1 Code Here
    end
end