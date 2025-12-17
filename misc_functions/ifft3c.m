function res = ifft3c(x)
%
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

res = sqrt(length(x(:)))*ifftshift(ifftn(fftshift(x)));