function Y = wavelet_per_frame(X, level)
%WAVELET_PER_FRAME  Apply wdenoise to each frame of [H W N]

X = single(X);
X = max(X,0);
[H,W,N] = size(X);

Y = zeros(H,W,N,'single');
parfor k = 1:N
    Y(:,:,k) = single(wdenoise(double(X(:,:,k)), level));
end

Y = max(Y,0);
end
