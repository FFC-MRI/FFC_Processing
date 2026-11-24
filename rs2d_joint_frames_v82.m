
function [Kout,info] = rs2d_joint_frames_v82(Kin)
% RS2D_JOINT_FRAMES_V82
%
% Background-guided model selection for residual coherent
% interference in 3D FCI k-space.
%
% Input dimensions:
% [kx ky kz evolutionTime evolutionField]
% [128 128 8 7 10]
%
% Requires existing rs2d_remove_coherent_tones.m.
%
% Experimental research code. Preserve original k-space
% and validate quantitative T1/R1 measurements separately.

assert(~isreal(Kin), 'Complex k-space required.');
assert(isequal(size(Kin),[128 128 8 7 10]), ...
    'Expected 128x128x8x7x10 FCI data.');
assert(all(isfinite(Kin(:))), ...
    'Input contains non-finite values.');

[nx,ny,nz,nt,nFields] = size(Kin);
nf = nt*nFields;
P = ny*nz;

allTones = [105.8204 64.5163 4.5315 17.8528 ...
            111.1902 11.1954 117.8495 124.5153];

targetX = 64.5163;

% Retain v8.1 limits for a controlled comparison.
nKnots = 16;
jointRank = 3;
ridge = 2;

minPEImprovement = 0.02;
minBackgroundImprovement = 0.45;

minGain = 0.05;
maxGain = 1.8;

maxGlobalChange = 0.06;
maxSliceChange = 0.12;
minStrength = 0.15;

modelNames = {'independent','joint35', ...
    'joint70','joint100','kz_adjusted'};
nModels = numel(modelNames);

%% 1. Existing first-stage correction

[Kbase,baseInfo] = rs2d_remove_coherent_tones( ...
    Kin, ...
    'SpatialDimensions',3, ...
    'CoilDimension',0, ...
    'PhaseOrder','centric-alternating', ...
    'KZFast',false, ...
    'ToneReadoutPositions',allTones, ...
    'RefineReadout',false, ...
    'MaxCandidates',18, ...
    'MaxTones',16, ...
    'EnvelopeKnots',24, ...
    'EnvelopeRidge',2, ...
    'MaxTotalAnatomyChange',0.91, ...
    'MaxPerSliceAnatomyChange',0.903, ...
    'HarmonicMode','auto');

Koriginal = reshape(double(Kin),[nx ny nz nf]);
K = reshape(double(Kbase),[nx ny nz nf]);

%% 2. Centric-alternating acquisition order

centre = floor(ny/2)+1;
order = centre;

for d = 1:ny
    if centre+d <= ny
        order(end+1) = centre+d; %#ok<AGROW>
    end
    if centre-d >= 1
        order(end+1) = centre-d; %#ok<AGROW>
    end
end

ranks = zeros(ny,1);
ranks(order) = 0:ny-1;

A = repmat(ranks,1,nz) + ...
    repmat((0:nz-1)*ny,ny,1);

pidx = double(A(:));

ky = (0:ny-1)';
tail = max(4,floor(0.30*ny));

outer = repmat(ky<tail | ky>=ny-tail,1,nz);
outer = outer(:);

train = outer & mod(pidx,4)~=3;
valid = outer & mod(pidx,4)==3;

%% 3. Smooth temporal envelope basis

knots = linspace(0,P-1,nKnots);
B = interp1(knots,eye(nKnots),pidx,'linear');

Bt = B(train,:);

D2 = diff(eye(nKnots),2,1);
G = Bt'*Bt;

scale = max(trace(G)/nKnots,eps);

R = G + ridge*scale*(D2'*D2) + ...
    1e-6*scale*eye(nKnots);

%% 4. Estimate residual phase and independent envelopes

n = (0:nx-1)'-floor(nx/2);

bx = exp(-2i*pi* ...
    (targetX-(floor(nx/2)+1))*n/nx);

Y = complex(zeros(P,nf));
theta = zeros(nf,1);
C = complex(zeros(nKnots,nf));

for f = 1:nf

    Y(:,f) = v82_project(K(:,:,:,f),bx);

    prior = v82_prior_theta(baseInfo,f,targetX);

    [theta(f),C(:,f)] = ...
        v82_fit_residual_theta( ...
        Y(:,f),prior,pidx,train,valid,B,Bt,R);
end

%% 5. Shared low-rank interference subspace

columnScale = sqrt(sum(abs(C).^2,1));

Cn = bsxfun(@rdivide,C,max(columnScale,eps));

[U,~,~] = svd(Cn,'econ');

r = min(jointRank,size(U,2));

U = U(:,1:r);
Cshared = U*(U'*C);

%% 6. Automatically construct anatomical support

support = v82_support(Koriginal);

if nnz(support)<150 || mean(support(:))>0.80
    error('v8.2: unreliable automatically generated support.');
end

background = convn(double(support), ...
    ones(5,5,1),'same')==0;

if nnz(background)<500
    background = ~support;
end

% Fixed spatial holdout pattern, shared by all candidate
% models in a frame.
[xx,yy,zz] = ndgrid(1:nx,1:ny,1:nz);

holdoutPattern = mod( ...
    floor((xx-1)/4) + ...
    floor((yy-1)/4) + ...
    floor((zz-1)/2),3)==0;

%% 7. Diagnostic structure

blank = struct( ...
    'frame',0, ...
    'model','', ...
    'status','', ...
    'readoutPosition',targetX, ...
    'theta',NaN, ...
    'phaseEncodeHoldout',NaN, ...
    'backgroundHoldout',NaN, ...
    'appliedHoldoutReduction',NaN, ...
    'gain',complex(NaN,NaN), ...
    'gainMagnitude',NaN, ...
    'gainPhaseDegrees',NaN, ...
    'fullAnatomyChange',NaN, ...
    'fullMaxSliceChange',NaN, ...
    'strength',0, ...
    'anatomyChange',NaN, ...
    'maxSliceChange',NaN, ...
    'backgroundReduction',NaN);

tests = repmat(blank,nf,1);
candidateTests = repmat(blank,nf*nModels,1);

accepted = zeros(nf,1);
strengths = zeros(nf,1);
selectedModel = strings(nf,1);

%% 8. Compare ALL models using held-out background

for f = 1:nf

    ph = exp(2i*pi*theta(f)*pidx);

    predictions = cell(1,nModels);
    peScores = -inf(1,nModels);

    blend = [0 0.35 0.70 1];

    % Four global temporal-envelope models.
    for j = 1:4

        coef = (1-blend(j))*C(:,f) + ...
            blend(j)*Cshared(:,f);

        predictions{j} = ph.*(B*coef);

        peScores(j) = v82_pe_score( ...
            Y(:,f),predictions{j},valid);
    end

    % Fifth candidate: kz-adjusted version of the
    % best global temporal-envelope model.
    [~,bestGlobal] = max(peScores(1:4));

    predictions{5} = v82_partition_model( ...
        predictions{bestGlobal},Y(:,f),train,ny,nz);

    peScores(5) = v82_pe_score( ...
        Y(:,f),predictions{5},valid);

    Iwork = v82_ifft3(K(:,:,:,f));
    Ioriginal = v82_ifft3(Koriginal(:,:,:,f));

    % Build all five candidate interference images.
    % All candidates will use exactly the same background ROI.
    Kmodels = cell(1,nModels);
    Imodels = cell(1,nModels);

    for j = 1:nModels

        Kmodels{j} = bsxfun(@times, ...
            reshape(bx,[nx 1 1]), ...
            reshape(predictions{j},[1 ny nz]));

        Imodels{j} = v82_ifft3(Kmodels{j});
    end

    roi = v82_common_background_roi( ...
        Imodels,background);

    % One background training/validation partition
    % for the entire five-model comparison.
    bgTrain = roi & ~holdoutPattern;
    bgValid = roi & holdoutPattern;

    candidates = repmat(blank,nModels,1);

    bestEligible = 0;
    bestScore = -inf;

    for j = 1:nModels

        rec = blank;
        rec.frame = f;
        rec.model = modelNames{j};
        rec.theta = theta(f);
        rec.phaseEncodeHoldout = peScores(j);

        if peScores(j)<minPEImprovement
            rec.status = 'poor_phase_encode_prediction';
            candidates(j) = rec;
            continue;
        end

        if nnz(bgTrain)<40 || nnz(bgValid)<24
            rec.status = 'insufficient_background';
            candidates(j) = rec;
            continue;
        end

        Imodel = Imodels{j};

        % Fit complex gain using background TRAINING pixels.
        m = Imodel(bgTrain);
        y = Iwork(bgTrain);

        power = sum(abs(m).^2);

        if power<eps
            rec.status = 'weak_model';
            candidates(j) = rec;
            continue;
        end

        alpha = sum(conj(m).*y)/power;

        rec.gain = alpha;
        rec.gainMagnitude = abs(alpha);
        rec.gainPhaseDegrees = rad2deg(angle(alpha));

        if abs(alpha)<minGain
            rec.status = 'absent_or_already_removed';
            candidates(j) = rec;
            continue;
        end

        if abs(alpha)>maxGain
            rec.status = 'gain_too_large';
            candidates(j) = rec;
            continue;
        end

        % Evaluate on common, independently withheld pixels.
        mv = Imodel(bgValid);
        yv = Iwork(bgValid);

        holdoutReduction = 1 - ...
            sum(abs(yv-alpha*mv).^2)/ ...
            max(sum(abs(yv).^2),eps);

        rec.backgroundHoldout = holdoutReduction;

        % Record full-strength anatomical effect.
        fullCandidate = Iwork-alpha*Imodel;

        [fullChange,fullSlice] = v82_anatomy_change( ...
            fullCandidate,Ioriginal,support);

        rec.fullAnatomyChange = fullChange;
        rec.fullMaxSliceChange = max(fullSlice);

        if holdoutReduction<minBackgroundImprovement
            rec.status = 'failed_background_holdout';
            candidates(j) = rec;
            continue;
        end

        % Determine anatomically permissible strength.
        strength = v82_permitted_strength( ...
            Iwork,Ioriginal,Imodel,alpha,support, ...
            maxGlobalChange,maxSliceChange);

        rec.strength = strength;

        if strength<minStrength
            rec.status = 'anatomy_change_limit';
            candidates(j) = rec;
            continue;
        end

        correctedImage = Iwork-strength*alpha*Imodel;

        [gc,sc] = v82_anatomy_change( ...
            correctedImage,Ioriginal,support);

        rec.anatomyChange = gc;
        rec.maxSliceChange = max(sc);

        % Background validation after the actually permitted
        % subtraction, not just the full-strength prediction.
        appliedReduction = 1 - ...
            sum(abs(yv-strength*alpha*mv).^2)/ ...
            max(sum(abs(yv).^2),eps);

        rec.appliedHoldoutReduction = appliedReduction;

        if appliedReduction<=0
            rec.status = 'no_applied_holdout_improvement';
            candidates(j) = rec;
            continue;
        end

        % Also require improvement across the entire
        % automatically identified background.
        e0 = sum(abs(Iwork(background)).^2);
        e1 = sum(abs(correctedImage(background)).^2);

        wholeReduction = 1-e1/max(e0,eps);

        rec.backgroundReduction = wholeReduction;

        if wholeReduction<=0
            rec.status = 'no_background_improvement';
            candidates(j) = rec;
            continue;
        end

        rec.status = 'eligible';
        candidates(j) = rec;

        % Model selection criterion:
        % best held-out suppression AFTER anatomy limits.
        if appliedReduction>bestScore
            bestScore = appliedReduction;
            bestEligible = j;
        end
    end

    % Save all five attempted models for later diagnosis.
    first = (f-1)*nModels+1;
    last = f*nModels;

    if bestEligible>0

        j = bestEligible;
        rec = candidates(j);
        rec.status = 'corrected';

        % Subtract selected complex model.
        K(:,:,:,f) = K(:,:,:,f) - ...
            rec.strength*rec.gain*Kmodels{j};

        accepted(f) = 1;
        strengths(f) = rec.strength;
        selectedModel(f) = string(rec.model);

        candidates(j) = rec;

        for k = 1:nModels
            if k~=j && strcmp(candidates(k).status,'eligible')
                candidates(k).status = 'eligible_not_selected';
            end
        end

        tests(f) = rec;

    else

        % No candidate could be safely accepted.
        % Save the most promising rejected model,
        % rather than leaving the frame unexplained.
        values = [candidates.backgroundHoldout];
        values(~isfinite(values)) = -inf;

        if all(values==-inf)
            values = [candidates.phaseEncodeHoldout];
            values(~isfinite(values)) = -inf;
        end

        [~,j] = max(values);

        tests(f) = candidates(j);
        selectedModel(f) = string(candidates(j).model);
    end

    candidateTests(first:last) = candidates;

end

%% 9. Output

Kout = cast(reshape(K,size(Kin)),'like',Kin);

info = struct();
info.version = '8.2-experimental';
info.status = 'completed';

info.targetReadoutPosition = targetX;
info.baseInfo = baseInfo;

% One final outcome per frame (70 records).
info.tests = tests;

% Every model attempted in every frame (350 records).
info.candidates = candidateTests;

info.acceptedPerFrame = accepted;
info.strengths = strengths;
info.selectedModel = selectedModel;
info.theta = theta;

info.totalAccepted = sum(accepted);
info.supportFraction = mean(support(:));

info.settings = struct( ...
    'jointRank',jointRank, ...
    'minPEImprovement',minPEImprovement, ...
    'minBackgroundImprovement',minBackgroundImprovement, ...
    'maxGlobalChange',maxGlobalChange, ...
    'maxSliceChange',maxSliceChange);

end

%% =========================================================
% LOCAL FUNCTIONS
% ==========================================================

function y = v82_project(k,b)
[nx,ny,nz] = size(k);
y = (conj(b).'*reshape(k,nx,ny*nz)).'/nx;
end

function I = v82_ifft3(k)
I = k;
for d=1:3
    I = ifftshift(I,d);
    I = ifft(I,[],d);
    I = fftshift(I,d);
end
end

function theta = v82_prior_theta(baseInfo,f,x)
theta = NaN;

if ~isfield(baseInfo,'toneTests') || ...
        isempty(baseInfo.toneTests)
    return;
end

tt = baseInfo.toneTests(:);

for j=1:numel(tt)
    if tt(j).frame==f && ...
            abs(tt(j).readoutPosition-x)<1 && ...
            isfinite(tt(j).phaseCyclesPerReadout)
        theta = tt(j).phaseCyclesPerReadout;
        return;
    end
end
end

function [theta,c] = v82_fit_residual_theta( ...
    y,prior,pidx,train,valid,B,Bt,R)

P = numel(y);
nfft = 16384;

temp = complex(zeros(P,1));
temp(pidx(train)+1) = y(train);

spectrum = abs(fft(temp,nfft));

seeds = zeros(1,3);

for k=1:3
    [~,idx] = max(spectrum);

    seeds(k) = (idx-1)/nfft;

    exclude = mod((idx-5:idx+5)-1,nfft)+1;
    spectrum(exclude) = 0;
end

if isfinite(prior)
    seeds(end+1) = prior;
end

bestScore = -inf;
theta = seeds(1);
c = complex(zeros(size(B,2),1));

for k=1:numel(seeds)

    times = pidx(train);
    signal = y(train);

    objective = @(th) -abs(sum( ...
        signal.*exp(-2i*pi*th*times)))^2;

    th = mod(fminbnd(objective, ...
        seeds(k)-2/nfft, ...
        seeds(k)+2/nfft),1);

    phase = exp(2i*pi*th*pidx);
    demodulated = y(train).*conj(phase(train));

    coef = R\(Bt'*demodulated);
    prediction = phase.*(B*coef);

    score = v82_pe_score(y,prediction,valid);

    if score>bestScore
        bestScore = score;
        theta = th;
        c = coef;
    end
end
end

function score = v82_pe_score(y,p,valid)
den = sum(abs(y(valid)).^2);
score = 1-sum(abs(y(valid)-p(valid)).^2)/ ...
    max(den,eps);
end

function adjusted = v82_partition_model( ...
    pred,y,train,ny,nz)

adjusted = pred;

for z=1:nz

    range = (z-1)*ny+(1:ny);
    indices = range(train(range));

    p = pred(indices);
    s = y(indices);

    power = sum(abs(p).^2);

    if power<eps
        continue;
    end

    % Regularise partition-specific coefficient towards 1.
    lambda = 0.5*power;

    gz = (sum(conj(p).*s)+lambda)/(power+lambda);

    delta = gz-1;

    if abs(delta)>0.5
        gz = 1+0.5*delta/abs(delta);
    end

    adjusted(range) = gz*pred(range);
end
end

function support = v82_support(K)
[nx,ny,nz,nf] = size(K);

sf = unique(round(linspace(1,nf,min(5,nf))));
a = zeros(nx,ny,nz);

for f=sf
    m = abs(v82_ifft3(K(:,:,:,f)));

    q = max(v82_quantile(m(:),0.995),eps);

    a = a+m/q/numel(sf);
end

a = movmedian(a,11,1);
a = movmedian(a,5,2);

a = v82_smooth(a,1,2.5);
a = v82_smooth(a,2,2.0);
a = v82_smooth(a,3,0.6);

support = a>0.26*v82_quantile(a(:),0.995);

if exist('bwconncomp','file')==2
    cc = bwconncomp(support,6);

    if cc.NumObjects>0
        [~,idx] = max(cellfun(@numel,cc.PixelIdxList));

        largest = false(size(support));
        largest(cc.PixelIdxList{idx}) = true;

        if nnz(largest)>=150
            support = largest;
        end
    end
end
end

function a = v82_smooth(a,dim,sigma)
radius = ceil(3*sigma);
u = -radius:radius;

g = exp(-0.5*(u/sigma).^2);
g = g/sum(g);

shape = [1 1 1];
shape(dim) = numel(g);

a = convn(a,reshape(g,shape),'same');
end

function q = v82_quantile(v,p)
v = sort(double(v(isfinite(v))));

if isempty(v)
    q = 0;
    return;
end

idx = max(1,min(numel(v), ...
    round(1+p*(numel(v)-1))));

q = v(idx);
end

function roi = v82_common_background_roi( ...
    modelImages,background)
% One ROI for all candidate models.
%
% Uses the maximum NORMALISED model magnitude, avoiding
% a comparison biased towards the strongest candidate.

roi = false(size(background));

if nnz(background)<100
    return;
end

combined = zeros(size(background));

for j=1:numel(modelImages)

    magnitude = abs(modelImages{j});

    rmsValue = sqrt(mean( ...
        magnitude(background).^2));

    if rmsValue<eps
        continue;
    end

    normalised = magnitude/rmsValue;

    combined = max(combined,normalised);
end

if ~any(combined(:)>0)
    return;
end

for p = [0.96 0.93 0.90]

    threshold = v82_quantile( ...
        combined(background),p);

    candidate = background & ...
        combined>=threshold & combined>0;

    if nnz(candidate)>=100
        roi = candidate;
        return;
    end
end
end

function strength = v82_permitted_strength( ...
    Iwork,Ioriginal,Imodel,alpha,support, ...
    maxGlobal,maxSlice)
% Search strongest allowable correction.
%
% A descending grid avoids assuming that a complex
% magnitude-change metric is strictly monotonic.

strength = 0;

for s = linspace(1,0,33)

    candidate = Iwork-s*alpha*Imodel;

    [gc,sc] = v82_anatomy_change( ...
        candidate,Ioriginal,support);

    if gc<=maxGlobal && max(sc)<=maxSlice
        strength = s;
        return;
    end
end
end

function [gc,sc] = v82_anatomy_change( ...
    candidate,reference,support)

m = abs(candidate);
r = abs(reference);

den = max(mean(r(support)),eps);

gc = mean(abs(m(support)-r(support)))/den;

nz = size(r,3);
sc = zeros(1,nz);

for z=1:nz
    mask = support(:,:,z);

    if ~any(mask(:))
        continue;
    end

    mm = m(:,:,z);
    rr = r(:,:,z);

    sliceDen = max(mean(rr(mask)),eps);

    sc(z) = mean(abs(mm(mask)-rr(mask)))/sliceDen;
end
end
