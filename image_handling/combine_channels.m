function [combined_images] = combine_channels(images,opts,obj)
%multicoil reconstruction
%   Detailed explanation goes here


images = images(:,:,:,:,:,opts);


tempdim = size(images);

%images = reshape(images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,size(images,6)]); % old debug code, enforce a specific dimensionality (with singletons if needed)

images_temp = reshape(images,tempdim(1),tempdim(2),[],size(images,6)); %reduce everything to 4 dimensions (X,Y, everything else, number of receivers).


if length(opts)>1

    n_channels = size(images_temp,4);
    if n_channels >1

        % Sum-of-squares combine
        combined_images(:,:,:) = rssq(images_temp(:,:,:,:),4);


        combined_images = reshape(combined_images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints]); %put everything back to its original dimensionality.
    end

else
    combined_images = images(:,:,:,:,:,1); %displays the individual channel data for debugging purposes.


end


end

