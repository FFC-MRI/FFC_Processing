function registered_images = register_images(imagestack)
%REGISTER_IMAGES Compatible interactive entry point; computation is separate.
    masks = ffc_registration_masks(imagestack);
    registered_images = ffc_register_images(imagestack,masks);
end
