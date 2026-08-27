function obj = apply_post_reconstruction_geometry(obj, geomT, includeMaps)
%APPLY_POST_RECONSTRUCTION_GEOMETRY Rotate/flip reconstructed arrays only.
%
% This function deliberately has no access to originalcomplexkspace or
% complexkspace. It operates only on reconstructed image products after the
% reconstruction core has returned.

    if nargin < 3
        includeMaps = false;
    end

    members = {'compleximage', 'complexcombined', 'magimage', 'phaseimage'};
    if includeMaps
        members = [members, {'T1Maps', 'R1Maps'}];
    end

    for memberIndex = 1:numel(members)
        name = members{memberIndex};
        if isobject(obj)
            exists = isprop(obj, name);
        else
            exists = isstruct(obj) && isfield(obj, name);
        end
        if ~exists || isempty(obj.(name))
            continue;
        end
        obj.(name) = apply_inplane_geometry(obj.(name), geomT);
    end
end
