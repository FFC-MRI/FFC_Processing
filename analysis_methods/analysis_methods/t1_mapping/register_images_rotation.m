function [ registered_images ] = register_images_rotation(imagestack)
%register_images Registers a stack of input images to a specified image
%   Detailed explanation goes here

dims = size(imagestack);
nbrow = size(imagestack,1);
nbcol = size(imagestack,2);

imagestack = reshape(imagestack,nbrow,nbcol,[]);
registered_images = zeros(size(imagestack));

for n=2:size(imagestack,3)
moving = imgaussfilt(abs(imagestack(:,:,n)));
fixed =imgaussfilt(abs(imagestack(:,:,1)));
 [optimizer,metric] = imregconfig('multimodal');   
 movingRegisteredDefault = imregister(moving,fixed,'rigid',optimizer,metric);   
 imshowpair(movingRegisteredDefault,fixed)
title('A: Default Registration')   
    
end

registered_images = imagestack;
% registered_images = reshape(registered_images,dims);

end

