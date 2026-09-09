function surfaces = cvi_build_surfaces(contours)
%CVI_BUILD_SURFACES Connect serial cvi42 SAX contours into open 3-D meshes.
%
% The mesh is only a geometric carrier for plane intersection; it is not
% intended as a diagnostic volumetric segmentation.  Rings are re-sampled
% by arc length and cyclically aligned between adjacent slices.

    types = {'LVEndo','LVEpi','RVEndo'};
    phases = unique([contours.PhaseIndex]);
    template = struct('Type','','PhaseIndex',NaN,'Vertices',[],'Faces',[], ...
        'SliceCount',0,'SliceIndices',[]);
    surfaces = repmat(template,0,1);
    nRing = 160;

    for p = phases
        for t = 1:numel(types)
            selected = contours(strcmp({contours.Type},types{t}) & ...
                [contours.PhaseIndex] == p);
            selected = selected(arrayfun(@(x) size(x.PointsPatientMm,1) >= 3,selected));
            if numel(selected) < 2, continue; end

            centreCells = arrayfun(@(x) mean(x.PointsPatientMm,1),selected, ...
                'UniformOutput',false);
            centres = vertcat(centreCells{:});

            % Prefer the explicit enhanced-DICOM stack index.  Surface
            % geometry is independent of whether the stack runs apex->base or
            % base->apex, but using the DICOM order prevents accidental
            % re-ordering by an arbitrary PCA sign and makes diagnostics
            % reproducible.  Fall back to centroid progression for older data.
            sliceIndices = [selected.SliceIndex];
            if numel(sliceIndices) == numel(selected) && ...
                    all(isfinite(sliceIndices)) && numel(unique(sliceIndices)) == numel(selected)
                [~,order] = sort(sliceIndices);
            else
                [~,~,v] = svd(centres - mean(centres,1),'econ');
                normal = v(:,1);
                location = centres * normal;
                [~,order] = sort(location);
            end
            selected = selected(order);

            rings = cell(numel(selected),1);
            for s = 1:numel(selected)
                ring = resampleClosed(selected(s).PointsPatientMm,nRing);
                if s > 1
                    ring = alignRing(rings{s-1},ring);
                end
                rings{s} = ring;
            end
            vertices = vertcat(rings{:});
            faces = zeros((numel(rings)-1)*nRing*2,3);
            row = 0;
            for s = 1:numel(rings)-1
                a0 = (s-1)*nRing;
                b0 = s*nRing;
                for k = 1:nRing
                    kn = mod(k,nRing)+1;
                    row = row+1; faces(row,:) = [a0+k,a0+kn,b0+k];
                    row = row+1; faces(row,:) = [a0+kn,b0+kn,b0+k];
                end
            end
            item = template;
            item.Type = types{t};
            item.PhaseIndex = p;
            item.Vertices = vertices;
            item.Faces = faces;
            item.SliceCount = numel(rings);
            item.SliceIndices = [selected.SliceIndex];
            surfaces(end+1,1) = item; %#ok<AGROW>
        end
    end
end

function ring = resampleClosed(points,n)
    points = double(points);
    if norm(points(1,:) - points(end,:)) > 1e-9
        points(end+1,:) = points(1,:);
    end
    segment = vecnorm(diff(points,1,1),2,2);
    keep = [true; segment > 1e-9];
    points = points(keep,:);
    if norm(points(1,:) - points(end,:)) > 1e-9
        points(end+1,:) = points(1,:);
    end
    distance = [0;cumsum(vecnorm(diff(points,1,1),2,2))];
    if distance(end) <= 0
        ring = repmat(points(1,:),n,1); return;
    end
    query = linspace(0,distance(end),n+1).'; query(end) = [];
    ring = [interp1(distance,points(:,1),query,'linear'), ...
            interp1(distance,points(:,2),query,'linear'), ...
            interp1(distance,points(:,3),query,'linear')];
end

function ring = alignRing(reference,ring)
    n = size(ring,1);
    best = inf; bestRing = ring;
    for reverse = 0:1
        candidate = ring;
        if reverse, candidate = flipud(candidate); end
        % Coarse-to-fine cyclic search keeps import time modest.
        step = max(1,floor(n/20));
        shifts = unique([0:step:n-1,n-1]);
        scores = zeros(size(shifts));
        for k = 1:numel(shifts)
            shifted = circshift(candidate,shifts(k),1);
            scores(k) = sum(vecnorm(reference-shifted,2,2).^2);
        end
        [~,idx] = min(scores);
        coarse = shifts(idx);
        fine = mod((coarse-step):(coarse+step),n);
        for sh = fine
            shifted = circshift(candidate,sh,1);
            score = sum(vecnorm(reference-shifted,2,2).^2);
            if score < best, best = score; bestRing = shifted; end
        end
    end
    ring = bestRing;
end
