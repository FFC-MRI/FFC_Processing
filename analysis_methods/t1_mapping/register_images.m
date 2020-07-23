function [ registered_images ] = register_images(imagestack)
%register_images Registers a stack of input images to a specified image
%   Detailed explanation goes here

dims = size(imagestack);
nbrow = size(imagestack,1);
nbcol = size(imagestack,2);
[u,v] = meshgrid(1:nbcol,1:nbrow);
imagestack = reshape(imagestack,nbrow,nbcol,[]);
registered_images = zeros(size(imagestack));
% 
for n=1:size(imagestack,3)
    [deltac,deltar] = fastreg((abs(imgaussfilt(abs(imagestack(:,:,n)),0.001))),(abs(imgaussfilt(abs(imagestack(:,:,1)),0.001))));
    exp_phase_ramp=exp(1i*2*pi*((u.*deltar)/nbcol+(v.*deltac)/nbrow));    %apply the translation in fourier space
    im_corrected = ifft2c(fft2c(imagestack(:,:,n)).*exp_phase_ramp);
    registered_images(:,:,n)=im_corrected;
end
%  f = waitbar(0,'Registering Images');
% for n=1:size(imagestack,3)
%    waitbar(n/size(imagestack,3),f)
%    FIXED = abs(imagestack(:,:,1));
%    MOVING = abs(imagestack(:,:,n));
%     fixedRefObj = imref2d(size(FIXED));
%     movingRefObj = imref2d(size(MOVING));
%     [~, tform] = registerImages_multimodal(abs(imagestack(:,:,n)),abs(imagestack(:,:,1)));
%     registered_images(:,:,n) = imwarp(MOVING, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
%    
% end
registered_images = reshape(registered_images,dims);
% delete(f)
end

