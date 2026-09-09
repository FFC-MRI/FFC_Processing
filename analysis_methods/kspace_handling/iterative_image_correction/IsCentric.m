function [outIm,outFft] = IsCentric(image)

% test if the image is a centred Fourier transform and if its FFT is as
% well. Centred FFT needs to use fftc and ifftc whereas the others can use
% fft and ifft.

image = squeeze(image);
imageTest = abs(image);
sze = size(image);
if length(sze) > 2
    error('The input must be a 2-D image.')
end
middle = round(sze/2);
width = round(sze/4);

% find if there is more power at the centre or at the edges
maskCentre = false(sze);
maskCentre(middle(1)-width(1) : middle(1)+width(1),...
           middle(2)-width(2) : middle(2)+width(2)) = true;
maskOuter = fftshift(maskCentre);

outIm = sum(abs(imageTest(maskCentre(:)))) > sum(abs(imageTest(maskOuter(:))));

% test the FFT of the image
imageFt = abs(ifft2(image));
outFft = sum(abs(imageFt(maskCentre(:)))) > sum(abs(imageFt(maskOuter(:))));