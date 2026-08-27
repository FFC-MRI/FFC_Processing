function res = ifft3c(x, n1, n2, n3)
%IFFT3C Orthonormal centred inverse FFT over dimensions 1, 2 and 3 only.
%
% Higher dimensions (time, field and receiver) are preserved. The previous
% ifftn implementation transformed every non-singleton dimension, including
% the receiver dimension in multichannel 3D data.

    if nargin < 2 || isempty(n1), n1 = size(x, 1); end
    if nargin < 3 || isempty(n2), n2 = size(x, 2); end
    if nargin < 4 || isempty(n3), n3 = size(x, 3); end

    x = ifftshift(x, 1);
    x = ifftshift(x, 2);
    x = ifftshift(x, 3);

    res = ifft(x, n1, 1);
    res = ifft(res, n2, 2);
    res = ifft(res, n3, 3);

    res = fftshift(res, 1);
    res = fftshift(res, 2);
    res = fftshift(res, 3);
    res = sqrt(n1*n2*n3) .* res;
end
