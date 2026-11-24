function test_ffc_central_band_rf()
% Check the source-derived order and model selection from an ambiguous alias.
rng(431);
nx = 96; ny = 128; nz = 3; nt = 4;
[x,y,z] = ndgrid((0:nx-1)-floor(nx/2),(0:ny-1)-floor(ny/2),0:nz-1);
clean = complex(zeros(nx,ny,nz,nt,'single'));
tone = clean;
ky = y(1,:,1);
sequence = [-6:63 -7:-1:-64];
rank = zeros(1,ny);
rank(sequence+65) = 0:ny-1;
read = exp(-2i*pi*x(:,1,1)*(34.27-1-floor(nx/2))/nx);
for t = 1:nt
    img = exp(-((x/26).^2+(y/22).^2+((z-1)/3).^2).^3) ...
        .*exp(1i*(.08*z+.01*x));
    img = single(img + .004*(randn(nx,ny,nz)+1i*randn(nx,ny,nz)));
    clean(:,:,:,t) = fftshift(fftn(ifftshift(img)));
    for k = 1:nz
        tone(:,:,k,t) = 8*exp(1i*(.31*t+.84*k))*read*exp(2i*pi*.1222*rank);
    end
end
raw = clean+tone;
[known,explicit] = rs2d_predict_coherent_rf(raw,'central-band',false);
[out,auto] = rs2d_predict_coherent_rf(raw,'centric2d',false);
assert(strcmp(auto.receivers.phaseOrder,'central-band'));
assert(isequal(out,known) && any(explicit.receivers.accepted(:)));
assert(norm(out(:)-clean(:))/norm(tone(:)) < .2);
[replay,~] = rs2d_predict_coherent_rf(raw,'centric2d',false,auto);
assert(isequal(out,replay));
fprintf('Central-band 3D order selection and replay passed.\n');
end
