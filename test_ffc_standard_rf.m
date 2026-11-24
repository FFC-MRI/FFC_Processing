function results = test_ffc_standard_rf()
%TEST_FFC_STANDARD_RF Exercise the standard stage and 2D linear safeguards.
% Run on the installed FFC Processing path, which supplies ifft2c et al.
rng(909,'twister');
nx=96; ny=90; nz=2; nt=3; nc=3; target=112;
[x,y] = ndgrid((0:nx-1)-floor(nx/2),(0:ny-1)-floor(ny/2));
im = exp(-((x/24).^2+((y+4)/23).^2).^3).*exp(1i*(.2+.01*x));
clean = complex(zeros(nx,ny,nz,nt,1,nc,'single'));
for c=1:nc
    for z=1:nz
        for t=1:nt
            I = single(im.*(1+.05*c+.08*z+.03*t)+.002*(randn(nx,ny)+1i*randn(nx,ny)));
            clean(:,:,z,t,1,c) = fftshift(fft2(ifftshift(I)));
        end
    end
end
r = (0:ny-1)-floor(ny/2);
bx = exp(-2i*pi*x(:,1)*(8.67-1-floor(nx/2))/nx);
by = exp(2i*pi*r*.2113);
tone = complex(zeros(size(clean),'single'));
for c=[1 2]
    for t=1:nt
        tone(:,:,1,t,1,c) = (6+c)*exp(1i*(.7*c+.5*t))*bx*by;
    end
end
% One image slice has the tone, one is clean, and receiver 3 is clean.
measured = clean+tone;
source = ffc_center_pad(measured,[nx target nz]);
scan = fixtureScan(measured,source);
sourceBefore = source;
cfg = ffc_rf_correction_config(scan);
assert(strcmp(cfg.phaseOrder,'linear') && cfg.is2D);
assert(isequal(cfg.phaseRange,[12 101]));
[corrected,info] = ffc_correct_rf(source,cfg);
assert(strcmp(info.status,'corrected'));
assert(isequal(source,sourceBefore) && isa(corrected,'single') && ~isreal(corrected));
assert(isequal(corrected(:,:,2,:,:,:),source(:,:,2,:,:,:)));
assert(isequal(corrected(:,:,:,:,:,3),source(:,:,:,:,:,3)));
padding = [1:11 102:target];
assert(isequal(corrected(:,padding,:,:,:,:),source(:,padding,:,:,:,:)));
actual = corrected(:,12:101,:,:,:,:);
results.syntheticLinearRemovedPct = 100*(1-norm(actual(:)-clean(:))^2/norm(tone(:))^2);
assert(results.syntheticLinearRemovedPct>95);
[replay,~] = ffc_correct_rf(source,cfg,info);
assert(isequal(replay,corrected));
% Check a real-valued spatially narrow MR structure well inside support.
% The phase evolution by itself is indistinguishable from a single RF tone.
object = clean(:,:,1,:,:,1);
voxel = complex(zeros(nx,ny,'single'));
voxel(48,37) = 7;
object = object + fftshift(fft2(ifftshift(voxel)));
[unchanged,weak] = rs2d_predict_coherent_rf(object,'linear',true);
assert(isequal(unchanged,object) && strcmp(weak.status,'no_correction'));
% The same stage runs in standard reconstruction, before coil combination.
built = ffc_reconstruct(scan,false);
assert(strcmp(built.rf_correction_info.status,'corrected'));
assert(isequal(built.originalcomplexkspace,sourceBefore));
assert(isequal(built.complexkspace,corrected));
warm = ffc_reconstruct(built,false);
assert(isequal(warm.magimage,built.magimage) && isequaln(warm.rf_correction_info,built.rf_correction_info));
% Changing the receiver selection and enabling whitening reuses the RF fit.
built.noise_whitening = true;
built.multichannel_recon = [2 1];
whitened = ffc_reconstruct(built,false);
assert(isequaln(whitened.rf_correction_info,built.rf_correction_info));
assert(isequal(whitened.whitening_info.selectedIndices,[2 1]));
assert(size(whitened.whitening_info.psi,1)==nc);
assert(all(isfinite(whitened.magimage(:))));
singleSlice = ffc_correct_rf(source(:,:,1,:,:,:),cfg);
assert(isequal(size(singleSlice),size(source(:,:,1,:,:,:))));
scan.param.TRANSFORM_PLUGIN = 'Spiral';
assert(strcmp(ffc_rf_correction_config(scan).status,'unsupported_trajectory'));
scan.param.TRANSFORM_PLUGIN = 'Sequential4D';
scan.is_svs = true;
assert(strcmp(ffc_rf_correction_config(scan).status,'spectroscopy'));
scan.is_svs = false;
scan.originalcomplexkspace = source(:,1,:,:,:,:);
assert(strcmp(ffc_rf_correction_config(scan).status,'insufficient_imaging_samples'));
fprintf('2D linear synthetic power removed %.2f%%. Standard reconstruction tests passed.\n', ...
    results.syntheticLinearRemovedPct);
end

function scan = fixtureScan(measured,source)
scan = struct('originalcomplexkspace',source,'rawdata',measured, ...
    'TwoDimensional',1,'param',struct('TRANSFORM_PLUGIN','Sequential4DCine'), ...
    'is_svs',false,'sequence','Cardiac Cine','bladeangles',[], ...
    'source_generation',1,'fft_size',size(source,1),'fft_size_3d',[], ...
    'multichannel_recon',[],'noise_whitening',false,'biascorrect','None', ...
    'rf_correction_info',[],'whitening_cache',[],'whitening_cache_key',[],'combination','RSS','fov',240,'fov_phase',240,'fov_3d',[], ...
    'window_function','None','window_size',1,'retain_coil_images',false, ...
    'retain_phase_images',false);
end
