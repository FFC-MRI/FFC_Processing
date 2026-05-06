function res = fft3c(x, n1, n2, n3)
%FFT3C Orthonormal centered 3D FFT over dims 1, 2 & 3 only.
%      Optional sizes n1,n2,n3 allow zero-pad/crop without padarray.
%
%      This preserves any higher dimensions, e.g. coils/time/field:
%          x = [X Y Z T F C]
%      transforms only X,Y,Z.

if nargin < 2 || isempty(n1), n1 = size(x,1); end
if nargin < 3 || isempty(n2), n2 = size(x,2); end
if nargin < 4 || isempty(n3), n3 = size(x,3); end

x = ifftshift(x, 1);
x = ifftshift(x, 2);
x = ifftshift(x, 3);

res = fft(x, n1, 1);
res = fft(res, n2, 2);
res = fft(res, n3, 3);

res = fftshift(res, 1);
res = fftshift(res, 2);
res = fftshift(res, 3);

res = (1/sqrt(n1*n2*n3)) .* res;
end