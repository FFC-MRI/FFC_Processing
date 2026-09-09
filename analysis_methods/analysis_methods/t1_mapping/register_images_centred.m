function [ registered_images ] = register_images_centred(imagestack)
%register_images Registers a stack of input images to a specified image
%   Detailed explanation goes here

dims = size(imagestack);
nbrow = size(imagestack,1);
nbcol = size(imagestack,2);
[u,v] = meshgrid(1:nbcol,1:nbrow);
imagestack = reshape(imagestack,nbrow,nbcol,[]);
registered_images = zeros(size(imagestack));

for n=1:size(imagestack,3)
    [deltac,deltar] = fastreg((abs(imagestack(:,:,n))),(abs(imagestack(:,:,1))));
    exp_phase_ramp=exp(1i*2*pi*((u.*deltar)/nbcol+(v.*deltac)/nbrow));    %apply the translation in fourier space
    im_corrected = ifft2(fft2(imagestack(:,:,n)).*exp_phase_ramp);
    registered_images(:,:,n)=im_corrected;
end

registered_images = reshape(registered_images,dims);

end

