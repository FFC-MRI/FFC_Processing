function [sliceIndex, view] = ffc_nearest_anatomical_slice(scan, planeName, patientPointMm, geomT)
%FFC_NEAREST_ANATOMICAL_SLICE Find the plane closest to a patient-space point.

    if nargin < 4 || isempty(geomT)
        geomT = eye(3);
    end
    view = ffc_anatomical_plane_geometry(scan, planeName, 1, geomT);
    relativePoint = double(patientPointMm(:)) - ...
        view.geometry.centerPatientMm;
    coordinate = dot(relativePoint, view.basisNormalPatient);
    sliceIndex = round(coordinate ./ view.spacingNormalMm + ...
        (view.nSlices + 1) ./ 2);
    sliceIndex = min(max(sliceIndex, 1), view.nSlices);
    view = ffc_anatomical_plane_geometry(scan, planeName, sliceIndex, geomT);
end
