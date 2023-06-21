classdef TotalVariationFilter < NoiseFilterTemplate
%Give the Name that will be compared to in the GUI and Checker
    properties
    processName = 'Total Variation';
    end


%Example Code
    methods
        function filtered_images = Filter(obj,tempimages,kernel)
            parfor n=1:size(tempimages,3)
                filtered_images(:,:,n) = TVL1denoise(tempimages(:,:,n),kernel,100);
            end
            disp('TOTAL');
            %disp(tempimages);
            %disp(kernel);
        end
        
    end
end