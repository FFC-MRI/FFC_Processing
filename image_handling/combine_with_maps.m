function combined_images = combine_with_maps(images, maps, opts, obj)
%COMBINE_WITH_MAPS  Roemer combine using precomputed sensitivity maps
%
% images: [X Y slices time field coils] complex
% maps:   [X Y slices 1 1 coils] complex (broadcasted over time/field)
% opts:   logical coil mask
%
% output: [X Y slices time field] complex

    images = images(:,:,:,:,:,opts);
    maps   = maps(:,:,:,:,:,opts);

    lambda = 0; % optional
    if isfield(obj,'coil_sens_lambda') && ~isempty(obj.coil_sens_lambda), lambda = obj.coil_sens_lambda; end

    num = sum(conj(maps) .* images, 6);
    den = sum(abs(maps).^2, 6) + lambda + eps(class(num));
    combined_images = num ./ den;
end
