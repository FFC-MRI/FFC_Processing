function [dB,a,b] = linear_correction(dB)

% dB contains a lilnear term a*t+b that needs to be removed
% since linear variations of dB do not affect the image (too
% much). So we find a and b from dB:
% B = epsilon + a*n + b
% mean(B) = m1 = a*n*(n+1)/2 + b*n
%mean(cumsum(B)) = m2 = a*(n+1)*(n+2)/6 + b*(n+1)/2

sdB = cumsum(dB);
m1 = mean(dB);
m2 = mean(sdB);
n = length(dB);
a = (m1 - 2*m2/(n+1))/((n+1)/2 - (n+2)/3);
b = m1 - a*(n+1)/2;
dB = dB(:) - (a*(1:n)'+b);