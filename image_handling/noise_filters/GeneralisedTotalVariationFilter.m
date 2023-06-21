classdef GeneralisedTotalVariationFilter < NoiseFilterTemplate

%Give the Name that will be compared to in the GUI and Checker
    properties
    processName = 'Generalised Total Variation';
    end


%Example Code
    methods
        function filtered_images = Filter(obj,tempimages,kernel)
            parfor n=1:size(tempimages,3)
                filtered_images(:,:,n) = imtgvsmooth(tempimages(:,:,n),kernel,kernel,200);
            end
            disp('GTOTAL');
        end
    end
end