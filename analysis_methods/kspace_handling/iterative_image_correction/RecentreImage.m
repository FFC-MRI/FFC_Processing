function imageCentred = RecentreImage(image)

% this function centres the image and its FFT to avoid problems at
% reconstruction.

sze = size(image);
imageSample = image(1:sze(1),1:sze(2),1); % only consider the first image
[outIm,outFft] = IsCentric(imageSample);
imageCentred = image;

if ~outIm
    imageCentred = fftshift(fftshift(image,1),2);
end

if ~outFft
    imageCentred = ifft(ifft(fftshift(fftshift(fft(fft(image,[],1),[],2),1),2),[],1),[],2);
end