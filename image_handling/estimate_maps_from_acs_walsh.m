function maps = estimate_maps_from_acs_walsh(kspace, opts, obj)
%ESTIMATE_MAPS_FROM_ACS_WALSH  Self-calibrated sensitivity maps from ACS k-space
%
% kspace: [X Y slices time field coils] complex (already padded to recon grid)
% opts: logical coil mask (dim 6 is coils)
% maps:  [X Y slices 1 1 coils] complex, normalised so sum|maps|^2 ~= 1

    % Select coils
    kspace = kspace(:,:,:,:,:,opts);

    % Force 6D shape: [X Y S T F C] even if some dims are singleton
    sz = size(kspace);
    sz(end+1:6) = 1;
    X = sz(1); Y = sz(2); S = sz(3); T = sz(4); F = sz(5); C = sz(6);
    kspace = reshape(kspace, [X Y S T F C]);

    % ---- tunables (good defaults for 128x128 in vivo) ----
    acs = 24;          % 24 or 32
    covSigma = 3;      % covariance smoothing sigma (pixels)
    maskFrac = 0.02;   % mask threshold = maskFrac * p95(RSS)
    acs = get_setting(obj, 'espirit_acs', acs);
    covSigma = get_setting(obj, 'walsh_cov_sigma', covSigma);
    maskFrac = get_setting(obj, 'coil_sens_mask_frac', maskFrac);

    maps = zeros(X,Y,S,1,1,C,'like',kspace);

    % ACS index ranges
    cx = floor(X/2) + 1;
    cy = floor(Y/2) + 1;
    hx = floor(acs/2);
    hy = floor(acs/2);
    xr = max(1, cx-hx) : min(X, cx+hx-1);
    yr = max(1, cy-hy) : min(Y, cy+hy-1);

    % Precompute denom smoother for masked smoothing
    filtSz = max(3, 2*ceil(3*covSigma)+1);

    for s = 1:S
        % ---- Build ACS-only k-space, averaged across (T,F) for stability ----
        % Kslice: [X Y T F C]
        Kslice = kspace(:,:,s,:,:,:);

        % Average across T and F *per k-space point* in the ACS region
        Kavg = zeros(X,Y,1,1,C,'like',kspace);          % keep 5D for ifft2c convenience
        for c = 1:C
            block = Kslice(xr,yr,:,:,c);               % [acs acs T F]
            Kavg(xr,yr,1,1,c) = mean(mean(block,4,'omitnan'),3,'omitnan');
        end

        % ---- Low-res coil images from ACS ----
        L = ifft2c(Kavg);                               % [X Y 1 1 C]
        L = reshape(L, [X Y C]);                        % [X Y C]

        % ---- Object mask from low-res RSS ----
        rss = sqrt(sum(abs(L).^2, 3));
        p95 = prctile(rss(:), 95);
        thr = max(eps(class(rss)), maskFrac * p95);
        mask = rss > thr;

        mask_f = imgaussfilt(single(mask), covSigma, 'FilterSize', filtSz);
        denom  = imgaussfilt(mask_f, covSigma, 'FilterSize', filtSz) + eps('single');

        % ---- Walsh covariance: R_ij(r) = smooth( L_i * conj(L_j) ) ----
        R = zeros(X,Y,C,C,'like',L);
        for i = 1:C
            Li = L(:,:,i);
            for j = 1:C
                Rij = Li .* conj(L(:,:,j));
                num = complex( ...
                    imgaussfilt(real(Rij).*mask_f, covSigma, 'FilterSize', filtSz), ...
                    imgaussfilt(imag(Rij).*mask_f, covSigma, 'FilterSize', filtSz) );
                R(:,:,i,j) = num ./ denom;
            end
        end

        % ---- Principal eigenvector per pixel ----
        Smap = zeros(X,Y,C,'like',L);

        for x = 1:X
            for y = 1:Y
                if mask(x,y)
                    M = squeeze(R(x,y,:,:));     % [C C]
                    M = (M + M')/2;              % enforce Hermitian
                    [V,D] = eig(M, 'vector');
                    [~,idx] = max(real(D));
                    v = V(:,idx);

                    % Fix arbitrary phase using the most stable receiver.
                    [~, refIndex] = max(abs(v));
                    v = v .* exp(-1i*angle(v(refIndex)));

                    Smap(x,y,:) = v;
                end
            end
        end

        % Normalise: sum|S|^2 ~= 1
        nrm = sqrt(sum(abs(Smap).^2, 3));
        Smap = Smap ./ (nrm + eps(class(nrm)));

        maps(:,:,s,1,1,:) = reshape(Smap, [X Y 1 1 1 C]);
    end
end


function value = get_setting(obj, name, defaultValue)
    value = defaultValue;
    if isstruct(obj) && isfield(obj, name) && ~isempty(obj.(name))
        value = obj.(name);
    elseif isobject(obj) && isprop(obj, name) && ~isempty(obj.(name))
        value = obj.(name);
    end
end
