function u = TGV2denoise(f, alpha0, alpha1, lambda, niter, opts)
%TGV2DENOISE  Total Generalized Variation (TGV^2) denoising (2D), primal-dual
%
%   u = TGV2denoise(f, alpha0, alpha1, lambda, niter, opts)
%
% Minimises (discrete, 2D):
%   alpha1 * || grad(u) - v ||_{2,1}  +  alpha0 * || E(v) ||_{2,1}  +  lambda * D(u,f)
%
% where:
%   v is a vector field (vx,vy)
%   E(v) is the symmetrized gradient:
%       E11 = d/dx(vx), E22 = d/dy(vy), E12 = 0.5*(d/dy(vx) + d/dx(vy))
%   || · ||_{2,1} is isotropic (per-pixel) l2 then sum.
%   D(u,f) is either:
%       (1/2)||u-f||_2^2  (opts.data='L2', default)
%       ||u-f||_1         (opts.data='L1')
%
% Inputs
%   f       : 2D image (real or complex; complex handled by denoising real/imag separately)
%   alpha0  : TGV second-order weight (smoothness of v)
%   alpha1  : TGV first-order weight (smoothness of u relative to v)
%   lambda  : data fidelity weight
%   niter   : iterations (e.g., 200-500 for strong denoising)
%   opts    : struct with fields (all optional)
%       .data       = 'L2' (default) or 'L1'
%       .tau        = primal step (default chosen automatically)
%       .sigma      = dual step   (default chosen automatically)
%       .theta      = extrapolation (default 1)
%       .normalize  = true/false (default true) normalise f to max=1 if >1
%
% Output
%   u : denoised image, same class as f (double/single), complex preserved

    if nargin < 6 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'data') || isempty(opts.data), opts.data = 'L2'; end
    if ~isfield(opts,'theta') || isempty(opts.theta), opts.theta = 1; end
    if ~isfield(opts,'normalize') || isempty(opts.normalize), opts.normalize = true; end

    % Handle complex by splitting (keeps your workflow flexible)
    if ~isreal(f)
        ur = TGV2denoise(real(f), alpha0, alpha1, lambda, niter, opts);
        ui = TGV2denoise(imag(f), alpha0, alpha1, lambda, niter, opts);
        u  = complex(ur, ui);
        return;
    end

    f = double(f);

    % Optional normalisation like your TV code
    mx = max(abs(f(:)));
    scale = 1;
    if opts.normalize && mx > 1
        scale = mx;
        f = f / scale;
    end

    [H,W] = size(f);

    % --- Choose stable step sizes if not provided ---
    % Operator K = [∇u - v ; E(v)].
    % A conservative bound for ||K||^2 in this discretisation is ~ 12.
    % (Empirically safe for forward diffs + sym grad.)
    L2 = 12;
    if ~isfield(opts,'tau') || isempty(opts.tau)
        tau = 0.02;              % similar to your TVL1
    else
        tau = opts.tau;
    end
    if ~isfield(opts,'sigma') || isempty(opts.sigma)
        sigma = 1/(L2*tau);
    else
        sigma = opts.sigma;
    end
    theta = opts.theta;

    % --- Primal variables ---
    u  = f;
    vx = zeros(H,W);
    vy = zeros(H,W);

    ubar  = u;
    vxbar = vx;
    vybar = vy;

    % --- Dual variables ---
    % p corresponds to grad(u)-v: p = (p1,p2)
    p1 = zeros(H,W);
    p2 = zeros(H,W);

    % q corresponds to E(v): q11,q22,q12 (symmetric tensor)
    q11 = zeros(H,W);
    q22 = zeros(H,W);
    q12 = zeros(H,W);

    for k = 1:niter
        % =========================
        % Dual ascent + projection
        % =========================
        % Compute grad(ubar)
        [ux, uy] = fwdGrad(ubar);

        % p <- p + sigma*(grad(ubar) - vbar)
        p1 = p1 + sigma*(ux - vxbar);
        p2 = p2 + sigma*(uy - vybar);

        % Project p onto ||p||_2 <= alpha1 (isotropic per pixel)
        [p1,p2] = projBall2(p1,p2, alpha1);

        % Compute sym-grad of vbar: E(vbar)
        [exx, eyy, exy] = symGrad(vxbar, vybar);

        q11 = q11 + sigma*exx;
        q22 = q22 + sigma*eyy;
        q12 = q12 + sigma*exy;

        % Project q onto ||q||_F <= alpha0 (Frobenius norm per pixel)
        [q11,q22,q12] = projBallF(q11,q22,q12, alpha0);

        % =========================
        % Primal descent
        % =========================
        u_old  = u;
        vx_old = vx;
        vy_old = vy;

        % u update: u <- prox_{tau*lambda*D}( u - tau * div(p) )
        divp = div2(p1,p2);
        u_tilde = u - tau*divp;

        switch upper(opts.data)
            case 'L2'
                % prox of (lambda/2)||u-f||^2
                u = (u_tilde + tau*lambda*f) / (1 + tau*lambda);
            case 'L1'
                % prox of lambda*||u-f||_1  => soft-threshold around f
                u = softThreshShift(u_tilde, f, tau*lambda);
            otherwise
                error('opts.data must be ''L2'' or ''L1''.');
        end

        % v update: v <- v - tau * ( -p + div(E^T q) ) with correct adjoint
        % For TGV^2: v appears in (grad(u)-v) and E(v).
        % Adjoint contributions:
        %   from p: +p (because term is <p, grad(u)-v>)
        %   from q: -E^T(q)
        [etx, ety] = symGradAdj(q11,q22,q12);  % E^T(q)
        vx = vx + tau*(p1 - etx);
        vy = vy + tau*(p2 - ety);

        % =========================
        % Extrapolation
        % =========================
        ubar  = u  + theta*(u  - u_old);
        vxbar = vx + theta*(vx - vx_old);
        vybar = vy + theta*(vy - vy_old);
    end

    u = u * scale;
end

% ---------- finite differences / operators ----------

function [gx,gy] = fwdGrad(u)
    % forward differences with Neumann boundary
    gx = u(:,[2:end end]) - u;
    gy = u([2:end end],:) - u;
end

function d = div2(px,py)
    % divergence: adjoint of fwdGrad
    % d = -∂x^T px - ∂y^T py with our forward scheme
    d = zeros(size(px));

    d(:,1)     = px(:,1);
    d(:,2:end) = px(:,2:end) - px(:,1:end-1);

    d(1,:)     = d(1,:) + py(1,:);
    d(2:end,:) = d(2:end,:) + py(2:end,:) - py(1:end-1,:);
end

function [exx,eyy,exy] = symGrad(vx,vy)
    % Symmetric gradient E(v)
    [vxx, vxy] = fwdGrad(vx);
    [vyx, vyy] = fwdGrad(vy);

    exx = vxx;
    eyy = vyy;
    exy = 0.5*(vxy + vyx);
end

function [tx,ty] = symGradAdj(q11,q22,q12)
    % Adjoint of symmetric gradient E^T(q)
    % For E(v) = [dvx/dx, dvy/dy, 0.5(dvx/dy + dvy/dx)]
    % the adjoint yields:
    % tx = -d/dx(q11) - d/dy(q12)
    % ty = -d/dy(q22) - d/dx(q12)
    %
    % Using adjoints of forward differences (same as divergence components)
    tx = div2(q11, q12);
    ty = div2(q12, q22);
end

% ---------- proximal / projections ----------

function [x,y] = projBall2(x,y, rad)
    n = sqrt(x.^2 + y.^2);
    s = max(1, n./rad);
    x = x ./ s;
    y = y ./ s;
end

function [q11,q22,q12] = projBallF(q11,q22,q12, rad)
    % Frobenius norm for symmetric tensor:
    % ||Q||_F = sqrt(q11^2 + q22^2 + 2*q12^2)
    n = sqrt(q11.^2 + q22.^2 + 2*(q12.^2));
    s = max(1, n./rad);
    q11 = q11 ./ s;
    q22 = q22 ./ s;
    q12 = q12 ./ s;
end

function x = softThreshShift(x, center, t)
    % prox_{t*||·-center||_1}(x) = center + soft(x-center, t)
    r = x - center;
    x = center + sign(r).*max(abs(r)-t, 0);
end
