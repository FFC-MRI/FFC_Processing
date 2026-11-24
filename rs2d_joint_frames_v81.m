
function [Kout,info] = rs2d_joint_frames_v81(Kin)
%RS2D_JOINT_FRAMES_V81 Joint residual interference correction.
%
% Experimental prototype for [128 128 8 7 10] FCI data.
% The existing v7.1 function provides the initial correction.
%
% v8.1 changes:
%   - Only refits residual tone at readout position 64.5163
%   - Fits phase progression from corrected (residual) k-space
%   - Uses original phase estimate only as an alternative seed
%   - Selects independent/joint envelopes on held-out PE lines
%   - Tests a regularised, kz-dependent complex envelope
%   - Validates background subtraction independently
%   - Records full-strength anatomy changes before limiting
%
% Original input is never modified in-place.
% No manually supplied anatomical mask is required.

assert(~isreal(Kin), 'Complex k-space required.');
assert(isequal(size(Kin),[128 128 8 7 10]), ...
    'This prototype requires [128 128 8 7 10].');

[nx,ny,nz,nt,nfField] = size(Kin);
nf = nt*nfField;
P = ny*nz;

allTones = [105.8204 64.5163 4.5315 17.8528 ...
            111.1902 11.1954 117.8495 124.5153];
targetX = 64.5163;

% Experimental settings
nKnots = 16;
jointRank = 3;
ridge = 2;
minPEImprovement = 0.02;
minBackgroundImprovement = 0.45;
maxGain = 1.8;
minGain = 0.05;
maxGlobalChange = 0.06;
maxSliceChange = 0.12;
minStrength = 0.15;

% =========================================================
% 1. Initial correction using the existing implementation
% =========================================================
[Kbase,baseInfo] = rs2d_remove_coherent_tones( ...
    Kin, ...
    'SpatialDimensions',3, ...
    'CoilDimension',0, ...
    'PhaseOrder','centric-alternating', ...
    'KZFast',false, ...
    'ToneReadoutPositions',allTones, ...
    'RefineReadout',false, ...
    'MaxCandidates',8, ...
    'MaxTones',6, ...
    'EnvelopeKnots',24, ...
    'EnvelopeRidge',2, ...
    'MaxTotalAnatomyChange',0.51, ...
    'MaxPerSliceAnatomyChange',0.53, ...
    'HarmonicMode','off');

Koriginal = reshape(double(Kin),[nx ny nz nf]);
K = reshape(double(Kbase),[nx ny nz nf]);

% =========================================================
% 2. Inferred RS2D acquisition order
% =========================================================
centre = floor(ny/2)+1;
order = centre;

for d=1:ny
    if centre+d<=ny
        order(end+1)=centre+d; %#ok<AGROW>
    end
    if centre-d>=1
        order(end+1)=centre-d; %#ok<AGROW>
    end
end

ranks = zeros(ny,1);
ranks(order)=0:ny-1;

A = repmat(ranks,1,nz) + ...
    repmat((0:nz-1)*ny,ny,1);
pidx = double(A(:));

ky = (0:ny-1)';
tail = max(4,floor(0.30*ny));

outer = repmat(ky<tail | ky>=ny-tail,1,nz);
outer = outer(:);

train = outer & mod(pidx,4)~=3;
valid = outer & mod(pidx,4)==3;

% =========================================================
% 3. Temporal envelope basis
% =========================================================
knots = linspace(0,P-1,nKnots);
B = interp1(knots,eye(nKnots),pidx,'linear');

Bt = B(train,:);
D2 = diff(eye(nKnots),2,1);

G = Bt'*Bt;
scale = max(trace(G)/nKnots,eps);

R = G + ridge*scale*(D2'*D2) + ...
    1e-6*scale*eye(nKnots);

% =========================================================
% 4. Project the residual tone in all 70 frames
% =========================================================
n = (0:nx-1)'-floor(nx/2);

bx = exp(-2i*pi* ...
    (targetX-(floor(nx/2)+1))*n/nx);

Y = complex(zeros(P,nf));
theta = zeros(nf,1);
C = complex(zeros(nKnots,nf));
independentPE = zeros(nf,1);

for f=1:nf
    Y(:,f) = v81_project(K(:,:,:,f),bx);

    prior = v81_prior_theta(baseInfo,f,targetX);

    [theta(f),C(:,f),independentPE(f)] = ...
        v81_fit_residual_theta( ...
            Y(:,f),prior,pidx,train,valid,B,Bt,R);
end

% =========================================================
% 5. Shared low-rank envelope subspace
% =========================================================
columnScale = sqrt(sum(abs(C).^2,1));

Cn = bsxfun(@rdivide,C,max(columnScale,eps));

[U,~,~] = svd(Cn,'econ');

r = min(jointRank,size(U,2));

% Project original complex coefficients into learned subspace.
U = U(:,1:r);
Cshared = U*(U'*C);

% =========================================================
% 6. Automatic anatomical/background support
% =========================================================
support = v81_support(Koriginal);

if nnz(support)<150 || mean(support(:))>0.80
    error('v8.1 could not construct a reliable brain support.');
end

background = convn(double(support), ...
    ones(5,5,1),'same')==0;

if nnz(background)<500
    background = ~support;
end

% =========================================================
% 7. Fit the residual model independently in each frame
% =========================================================
record = struct( ...
    'frame',0, ...
    'readoutPosition',targetX, ...
    'theta',NaN, ...
    'model','', ...
    'status','', ...
    'phaseEncodeHoldout',NaN, ...
    'backgroundHoldout',NaN, ...
    'gain',complex(NaN,NaN), ...
    'gainMagnitude',NaN, ...
    'gainPhaseDegrees',NaN, ...
    'fullAnatomyChange',NaN, ...
    'fullMaxSliceChange',NaN, ...
    'strength',0, ...
    'anatomyChange',NaN, ...
    'maxSliceChange',NaN, ...
    'backgroundReduction',NaN);

tests = repmat(record,nf,1);
accepted = zeros(nf,1);
strengths = zeros(nf,1);

[xx,yy,zz] = ndgrid(1:nx,1:ny,1:nz);
holdoutPattern = mod( ...
    floor((xx-1)/4)+floor((yy-1)/4)+ ...
    floor((zz-1)/2),3)==0;

for f=1:nf
    rec = record;
    rec.frame=f;
    rec.theta=theta(f);

    ph = exp(2i*pi*theta(f)*pidx);

    % Candidate envelope blends:
    % 1 = independent, 4 = entirely shared-subspace.
    blend = [0 0.35 0.70 1.0];
    predictions = cell(1,5);
    scores = -inf(1,5);

    for j=1:4
        coef = (1-blend(j))*C(:,f) + ...
            blend(j)*Cshared(:,f);

        predictions{j} = ph.*(B*coef);
        scores(j) = v81_pe_score( ...
            Y(:,f),predictions{j},valid);
    end

    % Optional kz-dependent fit, regularised towards unity.
    [~,bestGlobal] = max(scores(1:4));

    predictions{5} = v81_partition_model( ...
        predictions{bestGlobal},Y(:,f),train,ny,nz);

    scores(5) = v81_pe_score( ...
        Y(:,f),predictions{5},valid);

    [bestScore,best] = max(scores);
    rec.phaseEncodeHoldout=bestScore;

    modelNames = {'independent','joint35', ...
        'joint70','joint100','kz_adjusted'};

    rec.model=modelNames{best};

    if bestScore<minPEImprovement
        rec.status='poor_phase_encode_prediction';
        tests(f)=rec;
        continue;
    end

    pred = predictions{best};

    % Single complex k-space interference model.
    Kmodel = bsxfun(@times, ...
        reshape(bx,[nx 1 1]), ...
        reshape(pred,[1 ny nz]));

    Iwork = v81_ifft3(K(:,:,:,f));
    Ioriginal = v81_ifft3(Koriginal(:,:,:,f));
    Imodel = v81_ifft3(Kmodel);

    % Build background ROI from predicted interference,
    % not the observed image or brain intensity.
    roi = v81_background_roi(abs(Imodel),background);

    if nnz(roi)<100
        rec.status='insufficient_background';
        tests(f)=rec;
        continue;
    end

    bgTrain = roi & ~holdoutPattern;
    bgValid = roi & holdoutPattern;

    if nnz(bgTrain)<40 || nnz(bgValid)<24
        rec.status='insufficient_holdout';
        tests(f)=rec;
        continue;
    end

    modelTrain = Imodel(bgTrain);
    signalTrain = Iwork(bgTrain);

    modelPower = sum(abs(modelTrain).^2);

    if modelPower<eps
        rec.status='weak_residual_model';
        tests(f)=rec;
        continue;
    end

    alpha = sum(conj(modelTrain).*signalTrain) / ...
        modelPower;

    rec.gain=alpha;
    rec.gainMagnitude=abs(alpha);
    rec.gainPhaseDegrees=rad2deg(angle(alpha));

    if abs(alpha)<minGain
        rec.status='absent_or_already_removed';
        tests(f)=rec;
        continue;
    end

    if abs(alpha)>maxGain
        rec.status='gain_too_large';
        tests(f)=rec;
        continue;
    end

    mv = Imodel(bgValid);
    yv = Iwork(bgValid);

    bgFit = 1-sum(abs(yv-alpha*mv).^2)/ ...
        max(sum(abs(yv).^2),eps);

    rec.backgroundHoldout=bgFit;

    if bgFit<minBackgroundImprovement
        rec.status='failed_background_holdout';
        tests(f)=rec;
        continue;
    end

    % Record the full subtraction's effect even if rejected.
    fullCandidate = Iwork-alpha*Imodel;

    [fullChange,fullSlice] = v81_anatomy_change( ...
        fullCandidate,Ioriginal,support);

    rec.fullAnatomyChange=fullChange;
    rec.fullMaxSliceChange=max(fullSlice);

    % Find the largest permissible correction strength.
    % Grid search avoids assuming magnitude change is
    % strictly monotonic under complex subtraction.
    chosenStrength=0;

    for s=linspace(1,0,65)
        candidate = Iwork-s*alpha*Imodel;

        [gc,sc] = v81_anatomy_change( ...
            candidate,Ioriginal,support);

        if gc<=maxGlobalChange && ...
                max(sc)<=maxSliceChange
            chosenStrength=s;
            break;
        end
    end

    if chosenStrength<minStrength
        rec.status='anatomy_change_limit';
        tests(f)=rec;
        continue;
    end

    correctedImage = ...
        Iwork-chosenStrength*alpha*Imodel;

    % Independent global-background improvement check.
    e0 = sum(abs(Iwork(background)).^2);
    e1 = sum(abs(correctedImage(background)).^2);

    actualReduction = 1-e1/max(e0,eps);
    rec.backgroundReduction=actualReduction;

    if actualReduction<=0
        rec.status='no_background_improvement';
        tests(f)=rec;
        continue;
    end

    % Accept subtraction in complex k-space.
    K(:,:,:,f) = K(:,:,:,f) - ...
        chosenStrength*alpha*Kmodel;

    [gc,sc] = v81_anatomy_change( ...
        correctedImage,Ioriginal,support);

    rec.status='corrected';
    rec.strength=chosenStrength;
    rec.anatomyChange=gc;
    rec.maxSliceChange=max(sc);

    tests(f)=rec;
    accepted(f)=1;
    strengths(f)=chosenStrength;
end

% =========================================================
% 8. Results
% =========================================================
Kout = cast(reshape(K,size(Kin)),'like',Kin);

info = struct();
info.version='8.1-experimental';
info.status='completed';
info.tones=targetX;
info.baseInfo=baseInfo;
info.tests=tests;
info.acceptedPerFrame=accepted;
info.strengths=strengths;
info.theta=theta;
info.independentPE=independentPE;
info.supportMask=support;
info.totalAccepted=sum(accepted);

info.settings=struct( ...
    'jointRank',jointRank, ...
    'minPEImprovement',minPEImprovement, ...
    'minBackgroundImprovement',minBackgroundImprovement, ...
    'maxGlobalChange',maxGlobalChange, ...
    'maxSliceChange',maxSliceChange);

end

% =========================================================
% Helpers
% =========================================================

function y=v81_project(k,b)
    [nx,ny,nz]=size(k);
    y=(conj(b).'*reshape(k,nx,ny*nz)).'/nx;
end

function I=v81_ifft3(k)
    I=k;
    for d=1:3
        I=ifftshift(I,d);
        I=ifft(I,[],d);
        I=fftshift(I,d);
    end
end

function theta=v81_prior_theta(baseInfo,f,x)
    theta=NaN;

    if ~isfield(baseInfo,'toneTests') || ...
            isempty(baseInfo.toneTests)
        return;
    end

    tt=baseInfo.toneTests(:);

    for j=1:numel(tt)
        if tt(j).frame==f && ...
                abs(tt(j).readoutPosition-x)<1 && ...
                isfinite(tt(j).phaseCyclesPerReadout)
            theta=tt(j).phaseCyclesPerReadout;
            return;
        end
    end
end

function [theta,c,bestScore]=v81_fit_residual_theta( ...
    y,prior,pidx,train,valid,B,Bt,R)

    P=numel(y);
    nfft=16384;

    temp=complex(zeros(P,1));
    temp(pidx(train)+1)=y(train);

    spectrum=abs(fft(temp,nfft));

    % Candidate phase frequencies from the actual residual.
    seeds=zeros(1,3);

    for k=1:3
        [~,idx]=max(spectrum);
        seeds(k)=(idx-1)/nfft;

        exclude=mod((idx-5:idx+5)-1,nfft)+1;
        spectrum(exclude)=0;
    end

    % The original estimate is a prior, not a fixed answer.
    if isfinite(prior)
        seeds(end+1)=prior;
    end

    bestScore=-inf;
    theta=seeds(1);
    c=complex(zeros(size(B,2),1));

    for k=1:numel(seeds)
        t=pidx(train);
        yy=y(train);

        objective=@(th) -abs(sum( ...
            yy.*exp(-2i*pi*th*t)))^2;

        th=mod(fminbnd(objective, ...
            seeds(k)-2/nfft, ...
            seeds(k)+2/nfft),1);

        phase=exp(2i*pi*th*pidx);

        z=y(train).*conj(phase(train));
        coef=R\(Bt'*z);

        prediction=phase.*(B*coef);
        score=v81_pe_score(y,prediction,valid);

        if score>bestScore
            bestScore=score;
            theta=th;
            c=coef;
        end
    end
end

function score=v81_pe_score(y,p,valid)
    denominator=sum(abs(y(valid)).^2);
    score=1-sum(abs(y(valid)-p(valid)).^2)/ ...
        max(denominator,eps);
end

function adjusted=v81_partition_model( ...
    pred,y,train,ny,nz)

    adjusted=pred;

    for z=1:nz
        range=(z-1)*ny+(1:ny);
        indices=range(train(range));

        p=pred(indices);
        s=y(indices);

        power=sum(abs(p).^2);

        if power<eps
            continue;
        end

        % Regularise complex kz gain towards one.
        lambda=0.5*power;

        gz=(sum(conj(p).*s)+lambda)/ ...
            (power+lambda);

        % Prevent extreme partition-specific extrapolation.
        delta=gz-1;

        if abs(delta)>0.5
            gz=1+0.5*delta/abs(delta);
        end

        adjusted(range)=gz*pred(range);
    end
end

function support=v81_support(K)
    [nx,ny,nz,nf]=size(K);
    sf=unique(round(linspace(1,nf,min(5,nf))));

    a=zeros(nx,ny,nz);

    for f=sf
        m=abs(v81_ifft3(K(:,:,:,f)));
        q=max(v81_quantile(m(:),0.995),eps);
        a=a+m/q/numel(sf);
    end

    a=movmedian(a,11,1);
    a=movmedian(a,5,2);

    a=v81_smooth(a,1,2.5);
    a=v81_smooth(a,2,2.0);
    a=v81_smooth(a,3,0.6);

    support=a>0.26*v81_quantile(a(:),0.995);

    % Image Processing Toolbox is optional.
    if exist('bwconncomp','file')==2
        cc=bwconncomp(support,6);

        if cc.NumObjects>0
            [~,idx]=max(cellfun( ...
                @numel,cc.PixelIdxList));

            largest=false(size(support));
            largest(cc.PixelIdxList{idx})=true;

            if nnz(largest)>=150
                support=largest;
            end
        end
    end
end

function a=v81_smooth(a,dim,sigma)
    radius=ceil(3*sigma);
    u=-radius:radius;

    g=exp(-0.5*(u/sigma).^2);
    g=g/sum(g);

    shape=[1 1 1];
    shape(dim)=numel(g);

    a=convn(a,reshape(g,shape),'same');
end

function q=v81_quantile(v,p)
    v=sort(double(v(isfinite(v))));

    if isempty(v)
        q=0;
        return;
    end

    idx=max(1,min(numel(v), ...
        round(1+p*(numel(v)-1))));

    q=v(idx);
end

function roi=v81_background_roi(modelMag,background)
    roi=false(size(background));

    for p=[0.96 0.93 0.90]
        threshold=v81_quantile( ...
            modelMag(background),p);

        candidate=background & ...
            modelMag>=threshold & modelMag>0;

        if nnz(candidate)>=100
            roi=candidate;
            return;
        end
    end
end

function [gc,sc]=v81_anatomy_change( ...
    candidate,reference,support)

    m=abs(candidate);
    r=abs(reference);

    den=max(mean(r(support)),eps);

    gc=mean(abs( ...
        m(support)-r(support)))/den;

    nz=size(r,3);
    sc=zeros(1,nz);

    for z=1:nz
        mask=support(:,:,z);

        if ~any(mask(:))
            continue;
        end

        mm=m(:,:,z);
        rr=r(:,:,z);

        sliceDen=max(mean(rr(mask)),eps);

        sc(z)=mean(abs( ...
            mm(mask)-rr(mask)))/sliceDen;
    end
end
