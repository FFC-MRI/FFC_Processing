function obj = apply_reconstruction_geometry(obj, geomT, includeMaps)
%APPLY_RECONSTRUCTION_GEOMETRY Backwards-compatible wrapper.
% New code should call apply_post_reconstruction_geometry to make the strict
% reconstruction-then-display-transform ordering explicit.

    if nargin < 3
        includeMaps = false;
    end
    obj = apply_post_reconstruction_geometry(obj, geomT, includeMaps);
end
