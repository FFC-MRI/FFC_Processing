function Y = wavelet_per_frame(X, level, options)
%WAVELET_PER_FRAME  Apply wdenoise to each frame of [H W N]

if nargin < 3, options = struct(); end
X = single(X);
X = max(X,0);
[H,W,N] = size(X);

Y = zeros(H,W,N,'single');
if ffc_use_parallel(options,N)
    parfor k = 1:N
        Y(:,:,k) = single(wdenoise(double(X(:,:,k)), level));
    end
else
    for k = 1:N
        Y(:,:,k) = single(wdenoise(double(X(:,:,k)), level));
    end
end

Y = max(Y,0);
end
