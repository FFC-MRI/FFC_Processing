function [imageOut,maskOut] = autodespike(imageIn,threshold,extension)
%AUTODESPIKE remove transients and spikes from raw FCI k-spaces.
%   [imageOut,maskOut] = autodespike(imageIn)
%   [imageOut,maskOut] = autodespike(imageIn,threshold)
%   [imageOut,maskOut] = autodespike(imageIn,threshold,extension)
%    
%   imageIn: raw signal from the FCI imager (k-space, untouched).
%   threshold: optional, default is 50.
%   extension: defines how much more of the spike is removed. default are 6
%   voxels over the frequency encode direction.
%   imageOut: k-space without the spikes.
%   maskOut: mask generated to find the spikes.
%   This function works by analysing the derivative of the image along all
%   dimensions except frequency encoding. A threshold is determined to
%   clear out the spikes, the default value is 50 but this can be modified
%   as an option to the function.
% LB 02/06/23
% Licenced under GNU GPLv3 
% https://www.gnu.org/licenses/gpl-3.0.en.html

sze = size(imageIn);
imageOut = imageIn;
if length(sze)<8
    sze((length(sze)+1):8) = 1;                 % makes sure the input is considered as 8-dim
end
if nargin<2
    threshold = -50; %default value
else
    threshold = -abs(threshold);                % only negative values can be used.
end
if nargin<3
    extension = 6;
end
szeExt = sze + [0 2 0 0 0 2 2 0];               % TODO: does not consider multi-slices yet, should add the third dimension once we hace multi-slide FCI.
maskOut = zeros(size(imageIn));
for coil = 1:sze(8)                             % process each channel separately, as they may provide vastly different signals
    reducedImage = imageIn(:,:,1,1,1,:,:,coil);         % remove the navigator data
    normImage = reducedImage/std(reducedImage(:));      % normalising to facilitate histogram analysis    
    extImage = zeros(szeExt);
    extImage(:,2:end-1,1,1,1,2:end-1,2:end-1)=normImage;    
    varImage = diff(diff(diff(abs(extImage),2,2),2,6),2,7);     % process the derivative over all relevant dimensions
    mask = varImage<threshold;                                  % find the spikes
    el = strel('line',extension,90);            % extend the voxels found to the frequency encode direction to remove all the spike
    mask = imdilate(mask,el);
    imageOut(mask(:)) = 0;                      % remove the spikes
    maskOut(:,:,:,:,:,:,:,coil) = mask;         % and save the mask
end


