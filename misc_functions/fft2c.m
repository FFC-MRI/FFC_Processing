function res = fft2c(x, n1, n2)
%FFT2C Orthonormal centered 2D FFT over dims 1&2 only.
%      Optional sizes n1,n2 allow zero-pad/crop without padarray.

if nargin < 2 || isempty(n1), n1 = size(x,1); end
if nargin < 3 || isempty(n2), n2 = size(x,2); end

x = ifftshift(x, 1);
x = ifftshift(x, 2);

res = fft2(x, n1, n2);

res = fftshift(res, 1);
res = fftshift(res, 2);

res = (1/sqrt(n1*n2)) .* res;
end
