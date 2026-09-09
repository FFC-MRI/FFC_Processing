function obj = apply_fov_offset_phase_ramp(obj)
%APPLY_FOV_OFFSET_PHASE_RAMP Apply prescribed translations in k-space.
%
% This is a reconstruction operation, not an image rotation. It applies the
% scanner's physical FOV offsets as Fourier phase ramps to a fresh working
% copy of the immutable preprocessed k-space. No display rotation or flip is
% read or applied here.

    obj = correct_orientation(obj);
end
