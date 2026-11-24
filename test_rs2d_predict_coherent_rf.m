function results = test_rs2d_predict_coherent_rf()
%TEST_RS2D_PREDICT_COHERENT_RF Toolbox-free synthetic tests with known truth.
rng(4817,'twister');
n=96; ny=96; nz=6; nt=7; nf=3;
[x,y,z] = ndgrid((0:n-1)-floor(n/2),(0:ny-1)-floor(ny/2), ...
    ((0:nz-1)-(nz-1)/2)/(nz/2));
brain = exp(-((x/24).^2+((y-1)/26).^2+(.35*z).^2).^3);
brain = brain.*(1+.13*sin(x/3).*cos(y/4));
coreROI = abs(x)<=4 & abs(y-10)<=5 & abs(z)<.8;
whole = brain>.15;
phase = .22 + x*.006 + y*.003;
te = [.015 .032 .063 .112 .205 .37 .66];
t1 = [.07 .16 .31];
clean = complex(zeros(n,ny,nz,nt,nf,'single'));
for f=1:nf
    for t=1:nt
        clean(:,:,:,t,f) = brain.*exp(1i*phase).*(1-1.8*exp(-te(t)/t1(f)));
    end
end
noisy = clean + single(.004*(randn(size(clean))+1i*randn(size(clean))));
base = spatialTransform(noisy,false);
positions = [47.02, 74.43, 13.68];
omega = [.557, .354, .741];
rank = (0:ny-1).'-floor(ny/2);
rank = 2*abs(rank)-(rank>0);
B = exp(-2i*pi*((0:n-1).'-floor(n/2))*(positions-1-floor(n/2))/n);
toneK = complex(zeros(size(base),'single'));
for f=1:nf
    for t=1:nt
        projected = complex(zeros(3,ny,nz));
        for j=1:3
            speed = omega(j) + .0012*sin(.4*t+.8*f+(1:nz));
            intercept = [7 2.7 1.7];
            amplitude = intercept(j)*(.65+.5*rand(1,nz)).*exp(2i*pi*rand(1,nz));
            projected(j,:,:) = exp(2i*pi*rank*speed).*amplitude;
        end
        toneK(:,:,:,t,f) = reshape(B*reshape(projected,3,[]),n,ny,nz);
    end
end
raw = base + toneK;
rawCopy = raw;
[corrected,info] = rs2d_predict_coherent_rf(raw);
assert(strcmp(info.status,'corrected'));
assert(strcmp(info.receivers.phaseOrder,'centric-positive'));
assert(all(min(abs(info.receivers.readoutPixel-positions),[],1)<.15));
assert(isequal(size(corrected),size(raw)) && isa(corrected,'single') && ~isreal(corrected));
assert(isequal(raw,rawCopy));
[rerun,~] = rs2d_predict_coherent_rf(raw);
assert(isequal(rerun,corrected),'A fresh correction of raw must be deterministic.');
untouched = rs2d_predict_coherent_rf(base);
assert(isequal(untouched,base),'Clean reference should pass through exactly.');
[~, badInfo] = rs2d_predict_coherent_rf(complex(NaN(64,64,2),ones(64,64,2)));
assert(strcmp(badInfo.status,'nonfinite_input'));
[~, smallInfo] = rs2d_predict_coherent_rf(complex(zeros(16),ones(16)));
assert(strcmp(smallInfo.status,'unsupported_input'));
dual = cat(6,double(raw),2*double(raw));
[dualCorrected,dualInfo] = rs2d_predict_coherent_rf(dual,'centric-positive');
assert(isa(dualCorrected,'double') && numel(dualInfo.receivers)==2);
channelOne = dualCorrected(:,:,:,:,:,1);
channelTwo = dualCorrected(:,:,:,:,:,2);
assert(norm(channelTwo(:)-2*channelOne(:)) < 1e-9*norm(channelOne(:)));
% The same temporal waveform on a different physical ky traversal.
[~,positionToKy] = sort(rank);
ky = (0:ny-1).'-floor(ny/2);
[~,indices] = sort(2*abs(ky)-(ky<0));
rankOfKy = zeros(ny,1);
rankOfKy(indices) = 1:ny;
tone = toneK(:,positionToKy(rankOfKy),:,:,1);
cleanField = base(:,:,:,:,1);
[repair,alt] = rs2d_predict_coherent_rf(cleanField+tone,'centric-auto');
assert(strcmp(alt.receivers.phaseOrder,'centric-negative'));
assert(norm(repair(:)-cleanField(:)) < .35*norm(tone(:)));
[disabled,disabledInfo] = ffc_correct_rf(raw,struct('status','disabled','version',2));
assert(isequal(disabled,raw) && strcmp(disabledInfo.status,'disabled'));

restored = spatialTransform(corrected,true);
rawImage = spatialTransform(raw,true);
worstQuality = abs(meanRankError(restored,clean,whole));
rawQuality = abs(meanRankError(rawImage,clean,whole));
assert(max(worstQuality./rawQuality)<.35, 'Expected improvement in every synthetic frame.');
pct = (1-norm(corrected(:)-base(:))^2/norm(toneK(:))^2)*100;
assert(pct>92,'Expected >92%% interference removal on the stationary-tone phantom.');
results.knownTonePowerRemovedPct = pct;
results.rawTissueNrmsePct = 100*rawQuality;
results.correctedTissueNrmsePct = 100*worstQuality;
results.t1TrueMs = t1*1000;
results.t1CleanMs = roiFits(clean,coreROI,te)*1000;
results.t1RawMs = roiFits(rawImage,coreROI,te)*1000;
results.t1CorrectedMs = roiFits(restored,coreROI,te)*1000;
assert(max(abs(results.t1CorrectedMs-results.t1TrueMs))<5);
results.detectedReadoutPixels = info.receivers.readoutPixel;
results.acceptedFrameTones = sum(info.receivers.accepted,2);
% Deliberate extrapolation failure: inversion of tone phase in the unseen core.
changedTones = toneK;
changedTones(:,abs((0:ny-1)-floor(ny/2))<floor(.19*ny),:,:,:) = ...
    -changedTones(:,abs((0:ny-1)-floor(ny/2))<floor(.19*ny),:,:,:);
[stress, stressInfo] = rs2d_predict_coherent_rf(base+changedTones);
stressImage = spatialTransform(stress,true);
stressRaw = spatialTransform(base+changedTones,true);
results.centralPhaseJumpRawTissueNrmsePct = 100*meanRankError(stressRaw,clean,whole);
results.centralPhaseJumpCorrectedTissueNrmsePct = 100*meanRankError(stressImage,clean,whole);
results.centralPhaseJumpAcceptedFrameTones = sum(stressInfo.receivers.accepted,2);
disp(results);
end

function q = meanRankError(image, truth, roi)
a = reshape(image,[],size(image,4)*size(image,5));
b = reshape(truth,size(a));
q = sqrt(sum(abs(a(roi(:),:)-b(roi(:),:)).^2,1)./sum(abs(b(roi(:),:)).^2,1));
end

function result = roiFits(I,roi,time)
R = reshape(I,[],size(I,4),size(I,5));
result = zeros(1,size(I,5));
for f=1:size(I,5)
    s = mean(R(roi(:),:,f),1);
    ph = angle(sum(s(end-1:end)));
    signed = abs(s).*sign(real(s.*exp(-1i*ph)));
    R1 = fminbnd(@(r) profile(r,time(:),double(signed(:))),.2,50);
    result(f) = 1/R1;
end
end

function loss = profile(r,t,y)
D = [ones(size(t)) -exp(-t*r)];
a = D\y;
if any(a<0)
    a = max(a,0);
    if a(2)==0, a(1)=max(mean(y),0); end
    if a(1)==0, a(2)=max(D(:,2)'*y/(D(:,2)'*D(:,2)),0); end
end
loss = sum((y-D*a).^2);
end

function X = spatialTransform(X,inverse)
for d=1:3
    X = ifftshift(X,d);
    if inverse, X = ifft(X,[],d); else, X = fft(X,[],d); end
    X = fftshift(X,d);
end
end
