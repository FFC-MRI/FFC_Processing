function test_ffc_centric2d_rf()
%TEST_FFC_CENTRIC2D_RF Validate off-object corroboration and a central phase jump.
rng(7001,'twister');
nx=128; ny=110; nt=4;
[x,y] = ndgrid((0:nx-1)-floor(nx/2),(0:ny-1)-floor(ny/2));
signal = exp(-(y/17).^6).*exp(1i*.01*x);
clean = complex(zeros(nx,ny,1,nt,'single'));
for t=1:nt
    image = single(signal+.003*(randn(nx,ny)+1i*randn(nx,ny)));
    clean(:,:,1,t) = fftshift(fft2(ifftshift(image)));
end
ky = y(1,:);
rank = 2*abs(ky)-(ky>0);
bx = exp(-2i*pi*x(:,1)*(52.38-1-floor(nx/2))/nx);
tone = complex(zeros(size(clean),'single'));
for t=1:nt
    tone(:,:,1,t) = 7.5*exp(1i*.7*t)*bx*exp(2i*pi*rank*.364);
end
raw = clean+tone;
[fixed,info] = rs2d_predict_coherent_rf(raw,'centric2d',true);
assert(any(info.receivers.slice.centralSpatialConfirmed(:)), ...
    'An off-object component should corroborate the central prediction.');
ratio = norm(fixed(:)-clean(:))/norm(tone(:));
assert(ratio<.2);
[repeat,~] = rs2d_predict_coherent_rf(raw,'centric2d',true,info);
assert(isequal(fixed,repeat));

% The middle lines can have the opposite phase despite an excellent
% shoulder prediction. The spatial confirmation must refuse this fit.
changed = tone;
changed(:,abs(ky)<=floor(.10*ny),:,:) = -changed(:,abs(ky)<=floor(.10*ny),:,:);
raw = clean+changed;
[refused,info] = rs2d_predict_coherent_rf(raw,'centric2d',true);
assert(isequal(refused,raw) && strcmp(info.status,'no_correction'));
fprintf('Centric2D spatial confirmation passed. Injected-tone power reduction %.2f%%.\n',100*(1-ratio^2));
end
