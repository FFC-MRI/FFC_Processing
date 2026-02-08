function combined_images = combine_with_maps(images, maps, opts, obj)
%COMBINE_WITH_MAPS  Roemer combine using precomputed sensitivity maps
%
% images: [X Y slices time field coils] complex
% maps:   [X Y slices 1 1 coils] complex (broadcasts over time/field)
% opts:   logical mask (length nCoils) OR numeric indices into *current* coil dim
%
% Important convention:
%   In the refactored workflow, coil selection is applied once upstream and the
%   coil dimension is renumbered to 1..nSel. Therefore opts should refer to this
%   *local* coil index space (usually 1:nSel).

    nCoils = size(images, 6);

    % Parse selection in the local index space
    if isempty(opts) || (isscalar(opts) && isnumeric(opts) && opts == 1)
        sel = 1:nCoils;
    elseif islogical(opts)
        if numel(opts) ~= nCoils
            error('combine_with_maps:maskLengthMismatch', ...
                'Logical coil mask length (%d) must match coil dimension (%d).', numel(opts), nCoils);
        end
        sel = find(opts).';
    else
        sel = opts(:).';
    end

    if isempty(sel)
        error('combine_with_maps:emptySelection', 'No coils selected.');
    end
    if any(sel < 1) || any(sel > nCoils)
        error('combine_with_maps:coilIndexOutOfRange', ...
            'Selected coil indices [%s] exceed available coils 1..%d.', num2str(sel), nCoils);
    end

    % If exactly one coil is selected: return it directly (no sensitivity normalisation)
    if numel(sel) == 1
        combined_images = images(:,:,:,:,:,sel);
        return;
    end

    images = images(:,:,:,:,:,sel);
    maps   = maps(:,:,:,:,:,sel);

    lambda = 0;
    if isfield(obj,'coil_sens_lambda') && ~isempty(obj.coil_sens_lambda)
        lambda = obj.coil_sens_lambda;
    end

    num = sum(conj(maps) .* images, 6);
    den = sum(abs(maps).^2, 6) + lambda + eps(class(num));
    combined_images = num ./ den;
end
