function nativeMask = ffc_restore_native_plane_mask(displayMask, view)
%FFC_RESTORE_NATIVE_PLANE_MASK Undo display rotation and plane transposition.
%
% This helper is valid when the displayed plane has native dimension 3 as
% its normal, which is the geometry accepted by ImageReconCore.T1dispersion.

    if view.normalDimension ~= 3
        error('ffc_restore_native_plane_mask:NotAcquiredPlane', ...
            ['ROI dispersion currently requires the acquired image plane ', ...
             '(native dimension 3 normal).']);
    end
    inverseTransform = view.effectiveTransform.';
    untransformed = apply_inplane_geometry(logical(displayMask), ...
        inverseTransform);
    if view.rowDimension == 1 && view.columnDimension == 2
        nativeMask = untransformed;
    elseif view.rowDimension == 2 && view.columnDimension == 1
        nativeMask = untransformed.';
    else
        error('ffc_restore_native_plane_mask:InvalidDimensions', ...
            'The displayed plane does not map onto native dimensions 1 and 2.');
    end
end
