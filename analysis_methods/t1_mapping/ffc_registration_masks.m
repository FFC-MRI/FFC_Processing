function masks = ffc_registration_masks(images)
%FFC_REGISTRATION_MASKS Collect ROIs on the client before numerical work.
    masks = false(size(images,1),size(images,2),size(images,3));
    for slice = 1:size(images,3)
        fig = figure;
        cleanup = onCleanup(@() delete_if_valid(fig));
        imshow(im2single(abs(images(:,:,slice,1,1))),[]);
        title(sprintf('Slice %d: draw ROI, double-click to finish',slice));
        roi = drawfreehand();
        if ~isvalid(roi)
            error('ffc_registration_masks:Cancelled','Registration ROI selection was cancelled.');
        end
        mask = createMask(roi);
        if ~any(mask(:)), mask(:) = true; end
        masks(:,:,slice) = mask;
        clear cleanup
    end
end
function delete_if_valid(fig)
    if isgraphics(fig), delete(fig); end
end
