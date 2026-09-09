function output = apply_inplane_geometry(input, geomT)
%APPLY_INPLANE_GEOMETRY Apply an exact orthogonal transform to dimensions 1/2.
%
% geomT is the accumulated 3-by-3 image-plane transform used by the GUI.
% Only the eight lossless 90-degree rotation/reflection states are accepted.
% No interpolation, resampling or change to higher dimensions is performed.

    if isempty(input) || isempty(geomT)
        output = input;
        return;
    end

    if ~isnumeric(geomT) || size(geomT, 1) < 2 || size(geomT, 2) < 2
        error('apply_inplane_geometry:invalidTransform', ...
            'geomT must contain a numeric 2-by-2 in-plane transform.');
    end

    rawT = double(geomT(1:2, 1:2));
    T = round(rawT);
    if any(~isfinite(rawT(:))) || max(abs(rawT(:) - T(:))) > 1e-8 || ...
            any(sum(abs(T), 1) ~= 1) || any(sum(abs(T), 2) ~= 1)
        error('apply_inplane_geometry:notDiscreteOrthogonal', ...
            'The in-plane transform must be a signed permutation matrix.');
    end

    permutation = [2 1 3:ndims(input)];

    if isequal(T, [1 0; 0 1])
        output = input;
    elseif isequal(T, [0 1; -1 0])
        % 90 degrees clockwise in displayed image coordinates.
        output = flip(permute(input, permutation), 2);
    elseif isequal(T, [-1 0; 0 -1])
        output = flip(flip(input, 1), 2);
    elseif isequal(T, [0 -1; 1 0])
        % 90 degrees anticlockwise in displayed image coordinates.
        output = flip(permute(input, permutation), 1);
    elseif isequal(T, [-1 0; 0 1])
        output = flip(input, 2);
    elseif isequal(T, [1 0; 0 -1])
        output = flip(input, 1);
    elseif isequal(T, [0 1; 1 0])
        output = permute(input, permutation);
    elseif isequal(T, [0 -1; -1 0])
        output = flip(flip(permute(input, permutation), 1), 2);
    else
        error('apply_inplane_geometry:unsupportedTransform', ...
            'Unsupported in-plane transform.');
    end
end
