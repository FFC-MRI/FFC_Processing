
function [Kout,info] = rs2d_acquisition_phase_v9(Kin,varargin)
% RS2D_ACQUISITION_PHASE_V9
%
% Experimental cross-frame coherent-interference correction.
%
% Requires rs2d_joint_frames_v82.m.
% Input: complex [128 128 8 7 10] k-space.
%
% Uses v8.2 first, then attempts a further correction of
% the residual component at readout position 64.5163.
%
% New models:
%   independent       Individual envelope, shared phase slope
%   shared            Shared envelope, complex frame amplitude
%   phase_locked      Shared envelope, constrained frame phase
%   field_shared      Shared envelope across seven time points
%   mixed             Independent/shared-envelope combination
%
% Without supplied timestamps, the phase-locked model tests
% a HYPOTHESISED continuous stored-frame acquisition order.
% It is not proof of physical phase continuity.
%
% Research prototype: retain original data and independently
% validate quantitative T1/R1 results.

p = inputParser;
addParameter(p,'ReadoutTimesSeconds',[], ...
    @(x) isempty(x)||isnumeric(x));
addParameter(p,'InterferenceFrequencyHz',[], ...
    @(x) isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)));
parse(p,varargin{:});
opt = p.Results;

assert(~isreal(Kin),'Complex k-space required.');
assert(isequal(size(Kin),[128 128 8 7 10]), ...
    'Expected 128x128x8x7x10 data.');
assert(all(isfinite(Kin(:))), ...
    'Input contains non-finite values.');

[nx,ny,nz,nt,nFields] = size(Kin);
nf = nt*nFields;
P = ny*nz;

targetX = 64.5163;
nKnots = 16;
ridge = 2;

% Retain previous experimental acceptance limits.
minPE = 0.02;
minBackground = 0.45;
minGain = 0.05;
maxGain = 1.8;
maxGlobalChange = 0.06;
maxSliceChange = 0.12;
minStrength = 0.15;

% ---------------------------------------------------------
% 1. Existing v8.2 correction
% ---------------------------------------------------------
[Kbase,baseInfo] = rs2d_joint_frames_v82(Kin);

Koriginal = reshape(double(Kin),[nx ny nz nf]);
K = reshape(double(Kbase),[nx ny nz nf]);

% Use the same inferred acquisition index as the v7 engine.
assert(isfield(baseInfo,'baseInfo') && ...
       isfield(baseInfo.baseInfo,'acquisitionIndex'), ...
       'Missing acquisition-order information.');

A = double(baseInfo.baseInfo.acquisitionIndex(:));
assert(numel(A)==P,'Invalid acquisition-index dimensions.');

% ---------------------------------------------------------
% 2. Project residual onto target readout tone
% ---------------------------------------------------------
n = (0:nx-1)'-floor(nx/2);

bx = exp(-2i*pi* ...
    (targetX-(floor(nx/2)+1))*n/nx);

Y = complex(zeros(P,nf));

for f=1:nf
    Y(:,f) = v9_project(K(:,:,:,f),bx);
end

% Peripheral phase-encoding lines.
ky = (0:ny-1)';
tail = max(4,floor(0.30*ny));

outer = repmat(ky<tail | ky>=ny-tail,1,nz);
outer = outer(:);

train = outer & mod(A,4)~=3;
valid = outer & mod(A,4)==3;

% ---------------------------------------------------------
% 3. Estimate a shared temporal phase slope
% ---------------------------------------------------------
nfft = 16384;
temp = complex(zeros(P,nf));

temp(A(train)+1,:) = Y(train,:);

spectrum = fft(temp,nfft,1);
power = sum(abs(spectrum).^2,2);

[~,peak] = max(power);
theta = (peak-1)/nfft;

% ---------------------------------------------------------
% 4. Build temporal phase model
% ---------------------------------------------------------
times = opt.ReadoutTimesSeconds;

if ~isempty(times)
    assert(isequal(size(times),[P nf]), ...
        'ReadoutTimesSeconds must be [1024 70].');

    assert(~isempty(opt.InterferenceFrequencyHz), ...
        'Supply InterferenceFrequencyHz with timestamps.');

    assert(all(isfinite(times(:))), ...
        'Readout timestamps must be finite.');

    % Remove arbitrary global time origin.
    times = double(times)-double(times(1,1));

    oscillator = exp(2i*pi* ...
        opt.InterferenceFrequencyHz*times);

    timingMode = 'measured_timestamps';

else
    % Hypothesis only: assumes sequential stored frames.
    % Do not interpret this as verified acquisition timing.
    pseudoTime = bsxfun(@plus,A,(0:nf-1)*P);

    oscillator = exp(2i*pi*theta*pseudoTime);

    timingMode = 'stored_order_hypothesis';
end

% ---------------------------------------------------------
% 5. Temporal envelope representation
% ---------------------------------------------------------
knots = linspace(0,P-1,nKnots);

B = interp1(knots,eye(nKnots),A,'linear');
Bt = B(train,:);

D2 = diff(eye(nKnots),2,1);
G = Bt'*Bt;

scale = max(trace(G)/nKnots,eps);

R = G + ridge*scale*(D2'*D2) + ...
    1e-6*scale*eye(nKnots);

% Demodulate the residual from its hypothesised
% acquisition-time phase.
Ztrain = Y(train,:).*conj(oscillator(train,:));

% Individual frame envelopes.
Ci = R\(Bt'*Ztrain);

predIndependent = oscillator.*(B*Ci);

% ---------------------------------------------------------
% 6. Shared envelope fitted jointly across 70 frames
% ---------------------------------------------------------
[U,~,~] = svd(Ztrain,'econ');

cShared = R\(Bt'*U(:,1));
hShared = B*cShared;

den = max(sum(abs(hShared(train)).^2),eps);

g = (hShared(train)'*Ztrain)/den;

% Fix the global complex phase convention so the
% phase-locked model can use nonnegative real gains.
rotation = angle(sum(g.^2))/2;

if sum(real(g*exp(-1i*rotation)))<0
    rotation = rotation+pi;
end

hShared = hShared*exp(1i*rotation);
g = g*exp(-1i*rotation);

predShared = oscillator.*(hShared*g);

% This model enforces the hypothesised phase relationship.
% It must survive background validation to be accepted.
gLocked = max(real(g),0);
predLocked = oscillator.*(hShared*gLocked);

% ---------------------------------------------------------
% 7. Shared envelope within each field
% ---------------------------------------------------------
predField = complex(zeros(P,nf));

for field=1:nFields
    ids = (field-1)*nt+(1:nt);

    [Uf,~,~] = svd(Ztrain(:,ids),'econ');

    cf = R\(Bt'*Uf(:,1));
    hf = B*cf;

    denf = max(sum(abs(hf(train)).^2),eps);
    gf = (hf(train)'*Ztrain(:,ids))/denf;

    predField(:,ids) = ...
        oscillator(:,ids).*(hf*gf);
end

% A partially shared model.
predMixed = 0.5*predIndependent + ...
            0.5*predShared;

modelNames = {'independent','shared', ...
    'phase_locked','field_shared','mixed'};

nModels = numel(modelNames);

% ---------------------------------------------------------
% 8. Automatic anatomical/background support
% ---------------------------------------------------------
support = v9_support(Koriginal);

assert(nnz(support)>=150 && mean(support(:))<0.80, ...
    'Automatic brain support is unreliable.');

background = convn(double(support), ...
    ones(5,5,1),'same')==0;

if nnz(background)<500
    background = ~support;
end

[xx,yy,zz] = ndgrid(1:nx,1:ny,1:nz);

holdoutPattern = mod( ...
    floor((xx-1)/4)+floor((yy-1)/4)+ ...
    floor((zz-1)/2),3)==0;

% ---------------------------------------------------------
% 9. Diagnostic records
% ---------------------------------------------------------
blank = struct( ...
    'frame',0, ...
    'model','', ...
    'status','', ...
    'phaseEncodeHoldout',NaN, ...
    'backgroundHoldout',NaN, ...
    'appliedHoldout',NaN, ...
    'gain',complex(NaN,NaN), ...
    'gainMagnitude',NaN, ...
    'gainPhaseDegrees',NaN, ...
    'fullGlobalChange',NaN, ...
    'fullMaxSliceChange',NaN, ...
    'strength',0, ...
    'backgroundReduction',NaN);

tests = repmat(blank,nf,1);
candidates = repmat(blank,nf*nModels,1);

accepted = zeros(nf,1);
strengths = zeros(nf,1);

% ---------------------------------------------------------
% 10. Independently test candidate models in each frame
% ---------------------------------------------------------
for f=1:nf

    preds = { ...
        predIndependent(:,f), ...
        predShared(:,f), ...
        predLocked(:,f), ...
        predField(:,f), ...
        predMixed(:,f)};

    Iwork = v9_ifft3(K(:,:,:,f));
    Ioriginal = v9_ifft3(Koriginal(:,:,:,f));

    Kmodels = cell(1,nModels);
    Imodels = cell(1,nModels);

    peScore = zeros(1,nModels);

    for j=1:nModels
        peScore(j) = 1 - ...
            sum(abs(Y(valid,f)-preds{j}(valid)).^2)/ ...
            max(sum(abs(Y(valid,f)).^2),eps);

        Kmodels{j} = bsxfun(@times, ...
            reshape(bx,[nx 1 1]), ...
            reshape(preds{j},[1 ny nz]));

        Imodels{j} = v9_ifft3(Kmodels{j});
    end

    % Common ROI for fair model comparison.
    roi = v9_roi(Imodels,background);

    bgTrain = roi & ~holdoutPattern;
    bgValid = roi & holdoutPattern;

    bestCandidate = 0;
    bestReduction = -inf;
    frameTests = repmat(blank,nModels,1);

    for j=1:nModels

        rec = blank;
        rec.frame=f;
        rec.model=modelNames{j};
        rec.phaseEncodeHoldout=peScore(j);

        if peScore(j)<minPE
            rec.status='poor_phase_encoding_fit';
            frameTests(j)=rec;
            continue;
        end

        if nnz(bgTrain)<40 || nnz(bgValid)<24
            rec.status='insufficient_background';
            frameTests(j)=rec;
            continue;
        end

        modelImage = Imodels{j};

        m = modelImage(bgTrain);
        y = Iwork(bgTrain);

        modelPower = sum(abs(m).^2);

        if modelPower<eps
            rec.status='weak_model';
            frameTests(j)=rec;
            continue;
        end

        if strcmp(modelNames{j},'phase_locked')
            % No freely fitted complex phase rotation:
            % otherwise the continuity test is meaningless.
            alpha = max(0,real(sum(conj(m).*y))/modelPower);
        else
            alpha = sum(conj(m).*y)/modelPower;
        end

        rec.gain=alpha;
        rec.gainMagnitude=abs(alpha);
        rec.gainPhaseDegrees=rad2deg(angle(alpha));

        if abs(alpha)<minGain
            rec.status='absent_or_already_removed';
            frameTests(j)=rec;
            continue;
        end

        if abs(alpha)>maxGain
            rec.status='gain_too_large';
            frameTests(j)=rec;
            continue;
        end

        mv = modelImage(bgValid);
        yv = Iwork(bgValid);

        fullHoldout = 1 - ...
            sum(abs(yv-alpha*mv).^2)/ ...
            max(sum(abs(yv).^2),eps);

        rec.backgroundHoldout=fullHoldout;

        fullImage = Iwork-alpha*modelImage;

        [gc,sc] = v9_anatomy( ...
            fullImage,Ioriginal,support);

        rec.fullGlobalChange=gc;
        rec.fullMaxSliceChange=max(sc);

        if fullHoldout<minBackground
            rec.status='failed_background_holdout';
            frameTests(j)=rec;
            continue;
        end

        % Cumulative anatomical protection relative to
        % the originally acquired image.
        strength=0;

        for s=linspace(1,0,65)

            candidateImage = ...
                Iwork-s*alpha*modelImage;

            [change,sliceChange] = v9_anatomy( ...
                candidateImage,Ioriginal,support);

            if change<=maxGlobalChange && ...
                    max(sliceChange)<=maxSliceChange

                strength=s;
                break;
            end
        end

        if strength<minStrength
            rec.status='anatomy_change_limit';
            frameTests(j)=rec;
            continue;
        end

        candidateImage = ...
            Iwork-strength*alpha*modelImage;

        applied = 1 - ...
            sum(abs(yv-strength*alpha*mv).^2)/ ...
            max(sum(abs(yv).^2),eps);

        rec.appliedHoldout=applied;

        if applied<=0
            rec.status='no_applied_improvement';
            frameTests(j)=rec;
            continue;
        end

        e0 = sum(abs(Iwork(background)).^2);
        e1 = sum(abs(candidateImage(background)).^2);

        bgReduction = 1-e1/max(e0,eps);
        rec.backgroundReduction=bgReduction;

        if bgReduction<=0
            rec.status='no_background_improvement';
            frameTests(j)=rec;
            continue;
        end

        rec.status='eligible';
        rec.strength=strength;
        frameTests(j)=rec;

        if applied>bestReduction
            bestReduction=applied;
            bestCandidate=j;
        end
    end

    if bestCandidate>0
        j=bestCandidate;
        rec=frameTests(j);

        K(:,:,:,f) = K(:,:,:,f) - ...
            rec.strength*rec.gain*Kmodels{j};

        rec.status='corrected';
        frameTests(j)=rec;

        tests(f)=rec;
        accepted(f)=1;
        strengths(f)=rec.strength;

        for k=1:nModels
            if k~=j && strcmp(frameTests(k).status,'eligible')
                frameTests(k).status='eligible_not_selected';
            end
        end
    else
        % Select the most informative rejection record.
        values = [frameTests.backgroundHoldout];
        values(~isfinite(values))=-inf;

        if all(values==-inf)
            values=[frameTests.phaseEncodeHoldout];
            values(~isfinite(values))=-inf;
        end

        [~,j]=max(values);
        tests(f)=frameTests(j);
    end

    candidates((f-1)*nModels+(1:nModels))=frameTests;
end

% ---------------------------------------------------------
% 11. Output and diagnostics
% ---------------------------------------------------------
Kout = cast(reshape(K,size(Kin)),'like',Kin);

info = struct();
info.version='9.0-acquisition-phase';
info.status='completed';

info.previousVersion=baseInfo.version;
info.previousAccepted=baseInfo.totalAccepted;

info.timingMode=timingMode;
info.sharedTheta=theta;
info.sharedEnvelope=hShared;
info.complexFrameAmplitudes=g;

info.amplitudeByTimeField=reshape(abs(g),[nt nFields]);
info.phaseByTimeField=reshape(angle(g),[nt nFields]);

info.tests=tests;
info.candidates=candidates;
info.acceptedPerFrame=accepted;
info.strengths=strengths;
info.additionalAccepted=sum(accepted);

info.settings=struct( ...
    'minPE',minPE, ...
    'minBackground',minBackground, ...
    'maxGlobalChange',maxGlobalChange, ...
    'maxSliceChange',maxSliceChange);

end

% =========================================================
% Helper functions
% =========================================================

function y=v9_project(k,b)
[nx,ny,nz]=size(k);
y=(conj(b).'*reshape(k,nx,ny*nz)).'/nx;
end

function I=v9_ifft3(k)
I=k;
for d=1:3
    I=ifftshift(I,d);
    I=ifft(I,[],d);
    I=fftshift(I,d);
end
end

function support=v9_support(K)
[nx,ny,nz,nf]=size(K);
sf=unique(round(linspace(1,nf,min(5,nf))));

a=zeros(nx,ny,nz);

for f=sf
    m=abs(v9_ifft3(K(:,:,:,f)));
    a=a+m/max(v9_quantile(m(:),0.995),eps)/numel(sf);
end

a=movmedian(a,11,1);
a=movmedian(a,5,2);

a=v9_smooth(a,1,2.5);
a=v9_smooth(a,2,2.0);
a=v9_smooth(a,3,0.6);

support=a>0.26*v9_quantile(a(:),0.995);

if exist('bwconncomp','file')==2
    cc=bwconncomp(support,6);

    if cc.NumObjects>0
        [~,ix]=max(cellfun(@numel,cc.PixelIdxList));

        largest=false(size(support));
        largest(cc.PixelIdxList{ix})=true;

        if nnz(largest)>=150
            support=largest;
        end
    end
end
end

function a=v9_smooth(a,dim,sigma)
r=ceil(3*sigma);
u=-r:r;

g=exp(-0.5*(u/sigma).^2);
g=g/sum(g);

shape=[1 1 1];
shape(dim)=numel(g);

a=convn(a,reshape(g,shape),'same');
end

function q=v9_quantile(v,p)
v=sort(double(v(isfinite(v))));

if isempty(v)
    q=0;
    return;
end

idx=max(1,min(numel(v), ...
    round(1+p*(numel(v)-1))));

q=v(idx);
end

function roi=v9_roi(modelImages,background)
roi=false(size(background));

if nnz(background)<100
    return;
end

combined=zeros(size(background));

for j=1:numel(modelImages)
    m=abs(modelImages{j});
    rmsValue=sqrt(mean(m(background).^2));

    if rmsValue>eps
        combined=max(combined,m/rmsValue);
    end
end

for p=[0.96 0.93 0.90]
    th=v9_quantile(combined(background),p);

    candidate=background & combined>=th & combined>0;

    if nnz(candidate)>=100
        roi=candidate;
        return;
    end
end
end

function [gc,sc]=v9_anatomy(candidate,reference,support)
m=abs(candidate);
r=abs(reference);

gc=mean(abs(m(support)-r(support)))/ ...
    max(mean(r(support)),eps);

nz=size(r,3);
sc=zeros(1,nz);

for z=1:nz
    mask=support(:,:,z);

    if ~any(mask(:))
        continue;
    end

    mm=m(:,:,z);
    rr=r(:,:,z);

    sc(z)=mean(abs(mm(mask)-rr(mask)))/ ...
        max(mean(rr(mask)),eps);
end
end
