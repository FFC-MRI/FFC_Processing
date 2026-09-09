function registered_images = ffc_register_images(imagestack,masks,options)
%FFC_REGISTER_IMAGES Demons computation, with no graphics or pool startup.
    if nargin < 3, options = struct(); end
    dims = size(imagestack);
    nbslices = size(imagestack,3);
    imagestack = reshape(imagestack,size(imagestack,1),size(imagestack,2),nbslices,[]);
    if ~isequal(size(masks,1),dims(1)) || ~isequal(size(masks,2),dims(2)) || size(masks,3) ~= nbslices
        error('ffc_register_images:MaskSize','Registration masks must match the spatial dimensions.');
    end
registered_images = zeros(size(imagestack), 'like', imagestack);

% Demons parameters
N = [400 200 100];        % iterations per pyramid level
smoothField = 1.5;

for s = 1:nbslices

    FIXED = im2single(abs(imagestack(:,:,s,1)));

    mask = masks(:,:,s);

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

    if ffc_use_parallel(options,size(imagestack,4))
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
    else
        for n = 1:size(imagestack,4)

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

end

registered_images = reshape(registered_images, dims);


end
