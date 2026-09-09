function Y = llr_denoise_stack(X, opts)
%LLR_DENOISE_STACK  Locally Low-Rank denoising for a stack [X,Y,N]
%
% X: [H,W,N] (double/single)
% opts.patch : patch size (default 8)
% opts.step  : stride (default 3)
% opts.niter : passes (default 2)
% opts.tau   : shrink factor (default 2.2) => threshold = tau*sigma*sqrt(max(m,n))
% opts.sigma : noise std estimate in X (if empty, auto-estimate from highpass)
%
% Returns Y same size as X.

if nargin < 2, opts = struct(); end
if ~isfield(opts,'patch') || isempty(opts.patch), opts.patch = 8; end
if ~isfield(opts,'step')  || isempty(opts.step),  opts.step  = 3; end
if ~isfield(opts,'niter') || isempty(opts.niter), opts.niter = 2; end
if ~isfield(opts,'tau')   || isempty(opts.tau),   opts.tau   = 2.2; end
if ~isfield(opts,'sigma'), opts.sigma = []; end
if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = false; end

X = single(X);
[H,W,N] = size(X);

% Estimate sigma if not provided (robust high-pass MAD on a representative frame)
sigma = opts.sigma;
if isempty(sigma) || ~isfinite(sigma) || sigma <= 0
    k = max(1, round(N/2));
    sigma = estimate_sigma_hp(X(:,:,k));
end

p = opts.patch;
s = opts.step;

Y = X;
wgt = zeros(H,W,'single');
for i = 1:s:(H-p+1)
    for j = 1:s:(W-p+1)
        wgt(i:i+p-1,j:j+p-1) = wgt(i:i+p-1,j:j+p-1) + 1;
    end
end
covered = wgt(:) > 0;
thr = opts.tau * sigma * sqrt(max(p*p,N));

for it = 1:opts.niter
    acc = zeros(H,W,N,'single');

    for i = 1:s:(H-p+1)
        for j = 1:s:(W-p+1)
            patch = Y(i:i+p-1, j:j+p-1, :);          % [p,p,N]
            P = reshape(patch, p*p, N);              % [m,n]
            patchMean = mean(P,1,'omitnan');
            P = P - patchMean;             % remove per-frame mean (helps)

            % SVD shrinkage
            [U,Sv,V] = svd(P,'econ');
            sing = diag(Sv);

            % Threshold (Gavish-Donoho-ish scaling, simplified)
            sing_sh = max(sing - thr, 0);

            Psh = U * (sing_sh .* V');

            % add mean back
            Psh = Psh + patchMean;

            patch_out = reshape(Psh, p, p, N);

            acc(i:i+p-1, j:j+p-1, :) = acc(i:i+p-1, j:j+p-1, :) + patch_out;
        end
    end

    % Handle borders not covered by the loop by simple copy-through
    Ynew = reshape(Y, H*W, N);
    acc = reshape(acc, H*W, N);
    Ynew(covered,:) = acc(covered,:) ./ wgt(covered);
    Ynew = reshape(Ynew,H,W,N);

    % Light stabilisation: keep non-negative
    Ynew = max(Ynew, 0);

    Y = Ynew;

    if opts.verbose
        fprintf('LLR pass %d/%d, sigma~%.4g, thr~%.4g\n', it, opts.niter, sigma, opts.tau*sigma*sqrt(max(p*p,N)));
    end
end

end

function sigma = estimate_sigma_hp(im)
% Robust sigma estimate using a high-pass (Laplacian-like) and MAD
im = single(im);
hp = im - imgaussfilt(im, 1.0);
sigma = 1.4826 * median(abs(hp(:) - median(hp(:))));
sigma = max(sigma, 1e-6);
end
