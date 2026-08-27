function geometry = physical_image_display_geometry(fovRead, fovPhase, geomT, imageSize)
%PHYSICAL_IMAGE_DISPLAY_GEOMETRY Exact physical axes for a displayed image.
%
% fovRead and fovPhase are the full physical FOVs in the same units. The
% unrotated array has read samples in rows and phase samples in columns.
% geomT determines whether those axes have been exchanged by a 90-degree
% rotation or diagonal reflection. imageSize is the current displayed array
% size after that transform.

    fovRead = double(fovRead);
    fovPhase = double(fovPhase);
    if ~isscalar(fovRead) || ~isfinite(fovRead) || fovRead <= 0 || ...
            ~isscalar(fovPhase) || ~isfinite(fovPhase) || fovPhase <= 0
        error('physical_image_display_geometry:invalidFov', ...
            'Read and phase FOVs must be positive finite scalars.');
    end

    if numel(imageSize) < 2 || any(~isfinite(double(imageSize(1:2)))) || ...
            any(double(imageSize(1:2)) < 1)
        error('physical_image_display_geometry:invalidImageSize', ...
            'imageSize must contain positive row and column counts.');
    end
    nRows = double(imageSize(1));
    nColumns = double(imageSize(2));

    swapsAxes = false;
    if ~isempty(geomT)
        if ~isnumeric(geomT) || size(geomT, 1) < 2 || size(geomT, 2) < 2
            error('physical_image_display_geometry:invalidTransform', ...
                'geomT must contain a numeric 2-by-2 in-plane transform.');
        end
        inplaneT = round(double(geomT(1:2, 1:2)));
        swapsAxes = all(abs(inplaneT([1 4])) == 0) && ...
            all(abs(inplaneT([2 3])) == 1);
    end

    if swapsAxes
        fovX = fovRead;
        fovY = fovPhase;
    else
        fovX = fovPhase;
        fovY = fovRead;
    end

    pixelSizeX = fovX / nColumns;
    pixelSizeY = fovY / nRows;

    geometry = struct();
    geometry.swapsAxes = swapsAxes;
    geometry.fovX = fovX;
    geometry.fovY = fovY;
    geometry.pixelSizeXY = [pixelSizeX, pixelSizeY];
    geometry.xData = [-fovX/2 + pixelSizeX/2, ...
                       fovX/2 - pixelSizeX/2];
    geometry.yData = [-fovY/2 + pixelSizeY/2, ...
                       fovY/2 - pixelSizeY/2];
    geometry.xLim = [-fovX/2, fovX/2];
    geometry.yLim = [-fovY/2, fovY/2];
end
