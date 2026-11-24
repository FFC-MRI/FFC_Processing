
function [Kout, info] = rs2d_joint_frames_v8(Kin)
% RS2D_JOINT_FRAMES_V8
% Experimental joint multi-frame coherent interference removal.
%
% Dataset: [kx ky kz evolutionTime evolutionField]
%          [128 128 8 7 10], no coil dimension.
%
% Requires the existing rs2d_remove_coherent_tones.m.
% Does not require a manually supplied anatomical mask.
%
% IMPORTANT: Research prototype. Check image differences and
% T1/R1 bias before using corrected data quantitatively.

assert(~isreal(Kin), 'Expected complex k-space.');
assert(isequal(size(Kin), [128 128 8 7 10]), ...
    'v8 is configured for the existing brain dataset only.');

[nx,ny,nz,nt,nfield] = size(Kin);
nf = nt*nfield;
P = ny*nz;

tones = [105.8204 64.5163 4.5315 17.8528 ...
         111.1902 11.1954 117.8495 124.5153];

% Experimental safeguards.
maxGlobalChange = 0.06;
maxSliceChange = 0.12;
minBackgroundHoldout = 0.45;
maxComplexGain = 1.8;
minStrength = 0.15;
maxJointTonesPerFrame = 4;

% Joint-envelope settings.
nKnots = 16;
jointRank = 3;
jointBlend = 0.65;
envelopeRidge = 2;

% ---------------------------------------------------------
% Stage 1: Limited v7 temporal correction.
% Reserve more anatomical-change budget for the joint stage.
% The old harmonic stage is disabled.
% ---------------------------------------------------------
[Kbase, baseInfo] = rs2d_remove_coherent_tones( ...
    Kin, ...
    'SpatialDimensions', 3, ...
    'CoilDimension', 0, ...
    'PhaseOrder', 'centric-alternating', ...
    'KZFast', false, ...
    'ToneReadoutPositions', tones, ...
    'RefineReadout', false, ...
    'MaxCandidates', 8, ...
    'MaxTones', 6, ...
    'EnvelopeKnots', 24, ...
    'EnvelopeRidge', 2, ...
    'MaxTotalAnatomyChange', 0.01, ...
    'MaxPerSliceAnatomyChange', 0.03, ...
    'HarmonicMode', 'off');

Korig = reshape(double(Kin), [nx ny nz nf]);
K = reshape(double(Kbase), [nx ny nz nf]);

% ---------------------------------------------------------
% Acquisition-order mapping: centric alternating ky,
% complete ky sweep before advancing kz.
% ---------------------------------------------------------
c = floor(ny/2)+1;
ord = c;

for d = 1:ny
    if c+d <= ny
        ord(end+1) = c+d; %#ok<AGROW>
    end
    if c-d >= 1
        ord(end+1) = c-d; %#ok<AGROW>
    end
end

rank = zeros(ny,1);
rank(ord) = 0:ny-1;

A = repmat(rank,1,nz) + ...
    repmat((0:nz-1)*ny,ny,1);

pidx = double(A(:));

% Peripheral phase-encoding lines used for training.
tail = max(4,floor(0.30*ny));
ky = (0:ny-1)';
outer = ky<tail | ky>=ny-tail;
outer = repmat(outer,1,nz);
outer = outer(:);

train = outer & mod(pidx,4)~=3;
valid = outer & mod(pidx,4)==3;

% Shared temporal-envelope basis.
knots = linspace(0,P-1,nKnots);
B = interp1(knots,eye(nKnots),pidx,'linear');
Bt = B(train,:);

D2 = diff(eye(nKnots),2,1);
L = Bt'*Bt;
ridgeScale = max(trace(L)/nKnots,eps);

R = L + envelopeRidge*ridgeScale*(D2'*D2) + ...
    1e-6*ridgeScale*eye(nKnots);

% ---------------------------------------------------------
% Automatic anatomical support: only used internally.
% ---------------------------------------------------------
support = v8_support(Korig);
background = convn(double(support), ...
    ones(5,5,1),'same') == 0;

if nnz(background)<500
    background = ~support;
end

% Diagnostics for every tone in every frame.
blank = struct( ...
    'frame',0,'tone',0,'readoutPosition',NaN, ...
    'theta',NaN,'status','', ...
    'complexGain',complex(NaN,NaN), ...
    'phaseEncodeHoldout',NaN, ...
    'backgroundHoldout',NaN, ...
    'wholeBackgroundReduction',NaN, ...
    'strength',0,'anatomyChange',NaN, ...
    'maxSliceChange',NaN);

tests = repmat(blank,nf*numel(tones),1);
acceptedPerFrame = zeros(nf,1);
strengths = zeros(nf,numel(tones));
gain = complex(nan(nf,numel(tones)));

% ---------------------------------------------------------
% Stage 2: Fit each recurring tone jointly across 70 frames.
% ---------------------------------------------------------
for j = 1:numel(tones)

    x = tones(j);
    n = (0:nx-1)'-floor(nx/2);

    bx = exp(-2i*pi* ...
        (x-(floor(nx/2)+1))*n/nx);

    Y = complex(zeros(P,nf));
    theta = zeros(1,nf);

    % Estimate frame-specific phase progression.
    for f = 1:nf
        Y(:,f) = v8_project(K(:,:,:,f),bx);

        theta(f) = v8_get_theta( ...
            baseInfo,f,j,x, ...
            Korig(:,:,:,f),bx,pidx,train);
    end

    % Independent envelope estimates provide observations for
    % a shared low-rank envelope model.
    C = complex(zeros(nKnots,nf));

    for f = 1:nf
        ph = exp(2i*pi*theta(f)*pidx);
        demod = Y(train,f).*conj(ph(train));

        C(:,f) = R \ (Bt'*demod);
    end

    % Normalise columns so high-interference frames do not
    % completely dominate the shared subspace.
    columnScale = sqrt(sum(abs(C).^2,1));
    Cnormalized = bsxfun(@rdivide,C, ...
        max(columnScale,eps));

    [U,S,V] = svd(Cnormalized,'econ');
    r = min(jointRank,size(U,2));

    Cjoint = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    Cjoint = bsxfun(@times,Cjoint,columnScale);

    % Preserve an independent component to accommodate drift.
    Cfit = (1-jointBlend)*C + jointBlend*Cjoint;

    % -----------------------------------------------------
    % Apply separately calibrated complex corrections.
    % -----------------------------------------------------
    for f = 1:nf

        row = (j-1)*nf+f;
        rec = blank;
        rec.frame = f;
        rec.tone = j;
        rec.readoutPosition = x;
        rec.theta = theta(f);

        if acceptedPerFrame(f)>=maxJointTonesPerFrame
            rec.status = 'maximum_tones';
            tests(row) = rec;
            continue;
        end

        ph = exp(2i*pi*theta(f)*pidx);
        pred = ph.*(B*Cfit(:,f));

        % Independent phase-encoding holdout diagnostic.
        yval = Y(valid,f);
        pval = pred(valid);

        rec.phaseEncodeHoldout = 1 - ...
            sum(abs(yval-pval).^2)/ ...
            max(sum(abs(yval).^2),eps);

        Kmodel = bsxfun(@times, ...
            reshape(bx,[nx 1 1]), ...
            reshape(pred,[1 ny nz]));

        Iwork = v8_ifft3(K(:,:,:,f));
        Iref = v8_ifft3(Korig(:,:,:,f));
        Imodel = v8_ifft3(Kmodel);

        % Choose strong predicted ghost regions outside anatomy.
        roi = v8_background_roi(abs(Imodel),background);

        if nnz(roi)<80
            rec.status = 'insufficient_background';
            tests(row) = rec;
            continue;
        end

        [ix,iy,iz] = ndgrid(1:nx,1:ny,1:nz);

        holdout = roi & mod( ...
            floor((ix-1)/4)+floor((iy-1)/4)+ ...
            floor((iz-1)/2),3)==0;

        training = roi & ~holdout;

        if nnz(training)<40 || nnz(holdout)<24
            rec.status = 'insufficient_holdout';
            tests(row) = rec;
            continue;
        end

        m = Imodel(training);
        y = Iwork(training);

        alpha = sum(conj(m).*y)/ ...
            max(sum(abs(m).^2),eps);

        rec.complexGain = alpha;
        gain(f,j) = alpha;

        % No arbitrary phase limit: complex phase is estimated
        % from the background, then independently validated.
        if abs(alpha)>maxComplexGain || abs(alpha)<0.05
            rec.status = 'gain_out_of_bounds';
            tests(row) = rec;
            continue;
        end

        mv = Imodel(holdout);
        yv = Iwork(holdout);

        reduction = 1 - ...
            sum(abs(yv-alpha*mv).^2)/ ...
            max(sum(abs(yv).^2),eps);

        rec.backgroundHoldout = reduction;

        if reduction < minBackgroundHoldout
            rec.status = 'failed_background_holdout';
            tests(row) = rec;
            continue;
        end

        % Test subtraction against the ORIGINAL image.
        candidate = Iwork-alpha*Imodel;

        [globalChange,sliceChanges] = ...
            v8_anatomy_change(candidate,Iref,support);

        if globalChange<=maxGlobalChange && ...
                max(sliceChanges)<=maxSliceChange
            strength = 1;
        else
            lo = 0;
            hi = 1;

            for iteration = 1:12
                mid = (lo+hi)/2;
                candidate = Iwork-mid*alpha*Imodel;

                [gc,sc] = v8_anatomy_change( ...
                    candidate,Iref,support);

                if gc<=maxGlobalChange && ...
                        max(sc)<=maxSliceChange
                    lo = mid;
                else
                    hi = mid;
                end
            end

            strength = lo;
        end

        if strength<minStrength
            rec.status = 'anatomy_change_limit';
            tests(row) = rec;
            continue;
        end

        candidate = Iwork-strength*alpha*Imodel;

        % Require actual improvement across the background,
        % not just the selected calibration pixels.
        e0 = sum(abs(Iwork(background)).^2);
        e1 = sum(abs(candidate(background)).^2);

        wholeReduction = 1-e1/max(e0,eps);
        rec.wholeBackgroundReduction = wholeReduction;

        if wholeReduction<=0
            rec.status = 'no_background_improvement';
            tests(row) = rec;
            continue;
        end

        % Accepted subtraction in complex k-space.
        K(:,:,:,f) = K(:,:,:,f) - ...
            strength*alpha*Kmodel;

        [gc,sc] = v8_anatomy_change( ...
            candidate,Iref,support);

        rec.status = 'corrected';
        rec.strength = strength;
        rec.anatomyChange = gc;
        rec.maxSliceChange = max(sc);

        tests(row) = rec;

        acceptedPerFrame(f) = acceptedPerFrame(f)+1;
        strengths(f,j) = strength;
    end
end

Kout = cast(reshape(K,size(Kin)),'like',Kin);

info = struct();
info.version = '8.0-experimental';
info.status = 'completed';
info.tones = tones;
info.baseInfo = baseInfo;
info.tests = tests;
info.acceptedPerFrame = acceptedPerFrame;
info.strengths = strengths;
info.complexGains = gain;
info.totalAccepted = sum(acceptedPerFrame);
info.supportFraction = mean(support(:));
info.settings = struct( ...
    'jointRank',jointRank, ...
    'jointBlend',jointBlend, ...
    'nKnots',nKnots, ...
    'minBackgroundHoldout',minBackgroundHoldout, ...
    'maxGlobalChange',maxGlobalChange, ...
    'maxSliceChange',maxSliceChange);

end

% =========================================================
% Local helper functions
% =========================================================

function y = v8_project(k,b)
    [nx,ny,nz] = size(k);
    y = (conj(b).'*reshape(k,nx,ny*nz)).'/nx;
end

function I = v8_ifft3(k)
    I = k;
    for d=1:3
        I = ifftshift(I,d);
    end
    for d=1:3
        I = ifft(I,[],d);
    end
    for d=1:3
        I = fftshift(I,d);
    end
end

function theta = v8_get_theta( ...
    baseInfo,f,j,x,k,b,pidx,train)

    theta = NaN;

    if isfield(baseInfo,'toneTests') && ...
            ~isempty(baseInfo.toneTests)

        tt = baseInfo.toneTests(:);

        for q = 1:numel(tt)
            if tt(q).frame==f && tt(q).tone==j && ...
                    abs(tt(q).readoutPosition-x)<1 && ...
                    isfinite(tt(q).phaseCyclesPerReadout)

                theta = tt(q).phaseCyclesPerReadout;
                break;
            end
        end
    end

    if isfinite(theta)
        return;
    end

    y = v8_project(k,b);
    P = numel(y);
    nfft = 16384;

    temp = complex(zeros(P,1));
    temp(pidx(train)+1) = y(train);

    spectrum = fft(temp,nfft);
    [~,ind] = max(abs(spectrum));

    theta0 = (ind-1)/nfft;
    t = pidx(train);
    yy = y(train);

    objective = @(th) -abs(sum( ...
        yy.*exp(-2i*pi*th*t)))^2;

    theta = mod(fminbnd(objective, ...
        theta0-1/nfft,theta0+1/nfft),1);
end

function support = v8_support(K)
    [nx,ny,nz,nf] = size(K);

    sampleFrames = unique(round( ...
        linspace(1,nf,min(nf,5))));

    a = zeros(nx,ny,nz);

    for f = sampleFrames
        m = abs(v8_ifft3(K(:,:,:,f)));
        scale = max(v8_quantile(m(:),0.995),eps);

        a = a + m/scale/numel(sampleFrames);
    end

    a = movmedian(a,11,1);
    a = movmedian(a,5,2);

    a = v8_smooth(a,1,2.5);
    a = v8_smooth(a,2,2.0);
    a = v8_smooth(a,3,0.6);

    threshold = 0.26*v8_quantile(a(:),0.995);
    support = a>threshold;

    % Optional connectivity filtering if available.
    if exist('bwconncomp','file')==2
        cc = bwconncomp(support,6);

        if cc.NumObjects>0
            [~,largest] = max(cellfun( ...
                @numel,cc.PixelIdxList));

            cleaned = false(size(support));
            cleaned(cc.PixelIdxList{largest}) = true;

            if nnz(cleaned)>150
                support = cleaned;
            end
        end
    end
end

function a = v8_smooth(a,dim,sigma)
    r = ceil(3*sigma);
    u = -r:r;

    g = exp(-0.5*(u/sigma).^2);
    g = g/sum(g);

    shape = [1 1 1];
    shape(dim) = numel(g);

    a = convn(a,reshape(g,shape),'same');
end

function q = v8_quantile(v,p)
    v = sort(double(v(isfinite(v))));

    if isempty(v)
        q = 0;
        return;
    end

    index = max(1,min(numel(v), ...
        round(1+p*(numel(v)-1))));

    q = v(index);
end

function roi = v8_background_roi(modelMagnitude,background)
    roi = false(size(background));

    if nnz(background)<100
        return;
    end

    for pct = [0.96 0.93 0.90]
        threshold = v8_quantile( ...
            modelMagnitude(background),pct);

        candidate = background & ...
            modelMagnitude>=threshold & ...
            modelMagnitude>0;

        if nnz(candidate)>=80
            roi = candidate;
            return;
        end
    end
end

function [globalChange,sliceChange] = ...
    v8_anatomy_change(candidate,reference,support)

    m = abs(candidate);
    r = abs(reference);

    denominator = max(mean(r(support)),eps);

    globalChange = mean(abs( ...
        m(support)-r(support)))/denominator;

    nz = size(r,3);
    sliceChange = zeros(1,nz);

    for z = 1:nz
        mask = support(:,:,z);

        if ~any(mask(:))
            continue;
        end

        mm = m(:,:,z);
        rr = r(:,:,z);

        sliceDen = max(mean(rr(mask)),eps);

        sliceChange(z) = mean(abs( ...
            mm(mask)-rr(mask)))/sliceDen;
    end
end
