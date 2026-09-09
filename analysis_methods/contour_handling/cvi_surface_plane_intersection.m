function [x,y,segmentsPatient] = cvi_surface_plane_intersection(surface,view,registration,sourceCenter)
%CVI_SURFACE_PLANE_INTERSECTION Intersect registered CVI mesh with FCI plane.

    V = double(surface.Vertices);
    F = surface.Faces;
    sourceCenter = double(sourceCenter(:).');
    if isfield(registration,'targetCenterPatientMm') && ...
            numel(registration.targetCenterPatientMm) == 3
        targetCenter = double(registration.targetCenterPatientMm(:).');
    else
        targetCenter = double(view.geometry.centerPatientMm(:).');
    end
    T = double(registration.translationMm(:).');
    angles = deg2rad(double(registration.rotationDeg(:).'));
    cx=cos(angles(1)); sx=sin(angles(1));
    cy=cos(angles(2)); sy=sin(angles(2));
    cz=cos(angles(3)); sz=sin(angles(3));
    Rx=[1 0 0;0 cx -sx;0 sx cx];
    Ry=[cy 0 sy;0 1 0;-sy 0 cy];
    Rz=[cz -sz 0;sz cz 0;0 0 1];
    R = Rz*Ry*Rx;

    % Optional handedness correction. A true mirror reflection has
    % determinant -1 and therefore cannot be produced by any combination of
    % rigid rotations/translations.  Reflect about the registered source
    % centre using axes frozen to the selected FCI target view.  Keeping this
    % separate from R also means subsequent angle adjustments remain proper
    % rigid-body rotations.
    M = eye(3);
    if isfield(registration,'mirrorLeftRight') && logical(registration.mirrorLeftRight)
        if isfield(registration,'targetBasisXPatient') && ...
                numel(registration.targetBasisXPatient) == 3
            mirrorAxis = double(registration.targetBasisXPatient(:));
        else
            mirrorAxis = double(view.basisXPatient(:));
        end
        if norm(mirrorAxis) > 0
            mirrorAxis = mirrorAxis ./ norm(mirrorAxis);
            M = (eye(3) - 2.*(mirrorAxis*mirrorAxis.')) * M;
        end
    end
    if isfield(registration,'mirrorUpDown') && logical(registration.mirrorUpDown)
        if isfield(registration,'targetBasisYPatient') && ...
                numel(registration.targetBasisYPatient) == 3
            mirrorAxis = double(registration.targetBasisYPatient(:));
        else
            mirrorAxis = double(view.basisYPatient(:));
        end
        if norm(mirrorAxis) > 0
            mirrorAxis = mirrorAxis ./ norm(mirrorAxis);
            M = (eye(3) - 2.*(mirrorAxis*mirrorAxis.')) * M;
        end
    end

    % A separate in-plane control is easier to use for the common case in
    % which the clinical contours have the correct tilt but need rotating
    % within the FCI image.  The axis is fixed when the FCI registration
    % target is chosen, so changing the displayed scan does not change the
    % registration transform.
    if isfield(registration,'inPlaneRotationDeg') && ...
            isfinite(double(registration.inPlaneRotationDeg)) && ...
            abs(double(registration.inPlaneRotationDeg)) > eps
        if isfield(registration,'targetNormalPatient') && ...
                numel(registration.targetNormalPatient) == 3
            axisVector = double(registration.targetNormalPatient(:));
        else
            axisVector = double(view.basisNormalPatient(:));
        end
        if norm(axisVector) > 0
            axisVector = axisVector ./ norm(axisVector);
            a = deg2rad(double(registration.inPlaneRotationDeg));
            K = [0 -axisVector(3) axisVector(2); ...
                 axisVector(3) 0 -axisVector(1); ...
                 -axisVector(2) axisVector(1) 0];
            Rin = eye(3).*cos(a) + (1-cos(a)).*(axisVector*axisVector.') + sin(a).*K;
            R = Rin * R;
        end
    end
    V = (M * (V-sourceCenter).').';
    V = (R * V.').';
    V = V + targetCenter + T;

    p0 = double(view.slicePointPatientMm(:).');
    n = double(view.basisNormalPatient(:).'); n = n./norm(n);
    tol = 1e-5;
    x = []; y = []; segmentsPatient = zeros(0,3);
    displayCenter = double(view.geometry.centerPatientMm(:).');
    bx = double(view.basisXPatient(:).');
    by = double(view.basisYPatient(:).');

    for f = 1:size(F,1)
        tri = V(F(f,:),:);
        d = (tri-p0)*n.';
        points = zeros(0,3);
        edges = [1 2;2 3;3 1];
        for e = 1:3
            a=edges(e,1); b=edges(e,2);
            da=d(a); db=d(b);
            if abs(da) <= tol, points(end+1,:) = tri(a,:); end %#ok<AGROW>
            if da*db < -tol^2
                alpha = da/(da-db);
                points(end+1,:) = tri(a,:) + alpha*(tri(b,:)-tri(a,:)); %#ok<AGROW>
            elseif abs(db) <= tol
                points(end+1,:) = tri(b,:); %#ok<AGROW>
            end
        end
        if size(points,1) < 2, continue; end
        points = unique(round(points,6),'rows','stable');
        if size(points,1) < 2, continue; end
        % Choose the two most separated points if a triangle lies almost
        % exactly in the plane and contributes more than two candidates.
        if size(points,1) > 2
            best=[1 2]; bestD=-inf;
            for a=1:size(points,1)-1
                for b=a+1:size(points,1)
                    dd=sum((points(a,:)-points(b,:)).^2);
                    if dd>bestD, bestD=dd; best=[a b]; end
                end
            end
            points=points(best,:);
        end
        rel = points-displayCenter;
        xs = rel*bx.'; ys = rel*by.';
        x = [x;xs(1);xs(2);NaN]; %#ok<AGROW>
        y = [y;ys(1);ys(2);NaN]; %#ok<AGROW>
        segmentsPatient = [segmentsPatient;points;[NaN NaN NaN]]; %#ok<AGROW>
    end
end
