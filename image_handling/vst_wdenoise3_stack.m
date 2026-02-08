function Y = vst_wdenoise3_stack(X, opts)
%VST_WDENOISE3_STACK  Variance-stabilise magnitude MRI then 3D wavelet denoise.
% X: [H W N] magnitude stack (single/double)

if nargin<2, opts=struct(); end
if ~isfield(opts,'vst'), opts.vst=true; end
if ~isfield(opts,'method') || isempty(opts.method), opts.method='Bayes'; end
if ~isfield(opts,'denoise') || isempty(opts.denoise), opts.denoise='BlockJS'; end
if ~isfield(opts,'level'), opts.level=[]; end
if ~isfield(opts,'block'), opts.block=[]; end

X = single(X);
X = max(X,0);

% Optionally denoise in blocks along N to reduce motion/contrast smearing
if isempty(opts.block)
    blocks = {1:size(X,3)};
else
    bs = opts.block;
    n = size(X,3);
    blocks = arrayfun(@(k) k:min(k+bs-1,n), 1:bs:n, 'UniformOutput', false);
end

Y = zeros(size(X),'single');

for b = 1:numel(blocks)
    idx = blocks{b};
    V = X(:,:,idx);

    if opts.vst
        % Simple Rician-ish VST: Anscombe-like for magnitude.
        % Not perfect without sigma, but improves Gaussianity in low SNR.
        % If you know sigma, we can do a proper generalized Anscombe.
        Vt = sqrt(V + 3/8);
    else
        Vt = V;
    end

    % 3D wavelet denoise (treat N as third dimension)
    if isempty(opts.level)
        Vd = wdenoise3(Vt, 'DenoisingMethod', opts.denoise, 'ThresholdRule', opts.method);
    else
        Vd = wdenoise3(Vt, opts.level, 'DenoisingMethod', opts.denoise, 'ThresholdRule', opts.method);
    end

    if opts.vst
        % Approximate inverse of sqrt(x+3/8)
        Vout = max(Vd.^2 - 3/8, 0);
    else
        Vout = Vd;
    end

    Y(:,:,idx) = Vout;
end
end
