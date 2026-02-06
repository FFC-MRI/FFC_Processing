function [ registered_images ] = register_images(imagestack)
%register_images Registers a stack of input images to a specified image
%   Detailed explanation goes here

dims = size(imagestack);
nbrow = size(imagestack,1);
nbcol = size(imagestack,2);
nbslices = size(imagestack,3);
[u,v] = meshgrid(1:nbcol,1:nbrow);
 imagestack = reshape(imagestack,nbrow,nbcol,nbslices,[]);
registered_images = zeros(size(imagestack));
% 
% for f=1:nbslices
% for s=1:size(imagestack,4)
%     [deltac,deltar] = fastreg((abs(imgaussfilt(abs(imagestack(:,:,f,s)),0.001))),(abs(imgaussfilt(abs(imagestack(:,:,f,1)),0.001))));
%     exp_phase_ramp=exp(1i*2*pi*((u.*deltar)/nbcol+(v.*deltac)/nbrow));    %apply the translation in fourier space
%     im_corrected = ifft2c(fft2c(imagestack(:,:,f,s)).*exp_phase_ramp);
%     registered_images(:,:,f,s)=im_corrected;
% end
% end
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
% Assumes imagestack is at least [rows, cols, nbslices, nImgs]
% You already have: nbslices, dims

registered_images = zeros(size(imagestack), 'like', imagestack);

% Demons parameters
N = [400 200 100];        % iterations per pyramid level
smoothField = 1.5;

for s = 1:nbslices

    FIXED = im2single(abs(imagestack(:,:,s,1)));

    % ---- Hand-drawn mask (once per slice) ----
    figure; imshow(FIXED,[]);
    title(sprintf('Slice %d: draw ROI, double-click to finish', s));
    h = drawfreehand();
    mask = createMask(h);
    close;

    if ~any(mask(:))
        mask(:) = true;   % fallback: whole image
    end

    % ---- Crop to mask bounding box ----
    bb = regionprops(mask,'BoundingBox');
    bb = bb(1).BoundingBox;

    x1 = max(1, floor(bb(1)));
    y1 = max(1, floor(bb(2)));
    x2 = min(size(FIXED,2), ceil(bb(1)+bb(3)-1));
    y2 = min(size(FIXED,1), ceil(bb(2)+bb(4)-1));

    maskC  = mask(y1:y2, x1:x2);
    fixedC = FIXED(y1:y2, x1:x2);
    fixedC(~maskC) = 0;

    parfor n = 1:size(imagestack,4)

        MOVING = im2single(abs(imagestack(:,:,s,n)));
        movingC = MOVING(y1:y2, x1:x2);
        movingC(~maskC) = 0;

        % ---- Elastic registration ----
        [~, movingRegC] = imregdemons(movingC, fixedC, N, ...
            'AccumulatedFieldSmoothing', smoothField, ...
            'PyramidLevels', numel(N));

        % ---- Paste back ----
        regFull = MOVING;
        regFull(y1:y2, x1:x2) = movingRegC;

        registered_images(:,:,s,n) = cast(regFull, 'like', imagestack);
    end
end

registered_images = reshape(registered_images, dims);


% delete(f)
end

