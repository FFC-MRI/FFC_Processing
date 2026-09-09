function [I1, ph, A1, bkgd, jj] = iterative_images_correction_v9(A, thresh_bkgd, max_iterations, thresh_phase, known_bkgd)
%ITERATIVE_IMAGES_CORRECTION_V9  Fast ky-phase correction (linear ky ordering).
%
%   - avoids fmincon over Ny vars
%   - estimates ky phase curve using control points + smooth interpolation
%   - projects out constant + linear components (prevents global phase + shift)
%   - optional smoothness regularisation (handles drift + instability)

    I1 = [];
    ph = [];
    A1 = [];
    bkgd = [];
    jj = [];

    % Stack handling: reshape to [Nx Ny Npages] and process each page
    if ndims(A) > 2
        sze = size(A);
        A2 = reshape(A, size(A,1), size(A,2), []);
        nP = size(A2,3);

        I1 = zeros(size(A2), 'like', A2);
        A1 = zeros(size(A2), 'like', A2);
        ph = zeros(1, size(A2,2), nP, 'like', real(A2));

        parfor n = 1:nP
            [I1(:,:,n), ph(:,:,n), A1(:,:,n), ~, ~] = iterative_images_correction_v9( ...
                A2(:,:,n), thresh_bkgd, max_iterations, thresh_phase, known_bkgd);
        end

        I1 = reshape(I1, sze);
        A1 = reshape(A1, sze);
        bkgd = [];
        jj = [];
        return;
    end

    % ---- Single 2D k-space ----
    [Nx, Ny] = size(A);

    I1 = ifft2c(A);

    if nargin <= 4 || isempty(known_bkgd)
        bkgd = abs(I1) < (thresh_bkgd * max(abs(I1(:))));
    else
        bkgd = (known_bkgd == 1);
    end

    maskIdx = find(bkgd(:));

    ph = zeros(1, Ny, 'like', real(I1));
    A1 = A;

    % Tunables
    nCtrl_base       = max(8, min(32, round(Ny/8)));
    maxEvals_base    = 120;
    lambda_smooth    = 1e-2;
    updateMaskEvery  = 3;

    for jj = 1:max_iterations
        dph = optim_image_fast_reg(A1, Ny, maskIdx, nCtrl_base, maxEvals_base, lambda_smooth);

        ph = ph + dph;
        A1 = A .* exp(-1i * ph);  % implicit expansion along kx
        I1 = ifft2c(A1);

        if std(dph) < thresh_phase
            break;
        end

        if (nargin <= 4 || isempty(known_bkgd)) && mod(jj, updateMaskEvery) == 0
            bkgd = abs(I1) < (thresh_bkgd * max(abs(I1(:))));
            maskIdx = find(bkgd(:));
        end
    end
end


function dph = optim_image_fast_reg(A, Ny, maskIdx, nCtrl, maxEvals, lambda_smooth)
    ky = (1:Ny).';
    ctrlPos = round(linspace(1, Ny, nCtrl)).';

    % Remove constant + linear-in-ky components
    t = ky;
    G = [ones(Ny,1), t];
    P = eye(Ny) - G * ((G.'*G) \ G.');

    % Smoothness penalty (second difference)
    e = ones(Ny,1);
    D2 = spdiags([e -2*e e], 0:2, Ny-2, Ny);

    p0 = zeros(nCtrl,1);
    opts = optimset('Display','off', 'MaxFunEvals', maxEvals, 'MaxIter', maxEvals);

    pCtrl = fminsearch(@objfun, p0, opts);

    phi = interp1(ctrlPos, pCtrl, ky, 'pchip', 'extrap');
    phi = P * phi;
    dph = phi.';  % 1 x Ny

    function s = objfun(p)
        phi_ = interp1(ctrlPos, p, ky, 'pchip', 'extrap');
        phi_ = P * phi_;

        Ac = A .* exp(-1i * phi_.');
        Ic = abs(ifft2c(Ac));

        bkg_term = sum(Ic(maskIdx));
        smooth_term = lambda_smooth * sum((D2 * phi_).^2);

        s = bkg_term + smooth_term;
    end
end
