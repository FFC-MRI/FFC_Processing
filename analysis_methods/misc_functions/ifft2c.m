function res = ifft2c(x, n1, n2)
%IFFT2C Orthonormal centered 2D IFFT over dims 1&2 only.
%       Optional sizes n1,n2 allow zero-pad/crop without padarray.

if nargin < 2 || isempty(n1), n1 = size(x,1); end
if nargin < 3 || isempty(n2), n2 = size(x,2); end

% Convert centred k-space to MATLAB's unshifted FFT ordering. Using
% ifftshift here is important for odd matrix sizes; fftshift is only
% interchangeable with it for even dimensions.
x = ifftshift(x, 1);
x = ifftshift(x, 2);

res = ifft2(x, n1, n2);

res = fftshift(res, 1);
res = fftshift(res, 2);

% Orthonormal scaling for the 2D plane only
res = sqrt(n1*n2) .* res;
end
