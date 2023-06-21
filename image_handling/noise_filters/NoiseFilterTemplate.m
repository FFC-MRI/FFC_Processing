classdef (Abstract) NoiseFilterTemplate

%Give the Name that will be compared to in the GUI and Checker
    properties (Abstract)
    processName;
    end


%Example Code
    methods (Abstract)
        filtered_images = Filter(tempimages,kernel)
            %Insert Desired Filter Code Here
            %parfor n=1:size(tempimages,3)
            %filtered_images(:,:,n) = imtgvsmooth(tempimages(:,:,n),kernel,kernel,200)
%             disp('HERE');
%         end
    end
end