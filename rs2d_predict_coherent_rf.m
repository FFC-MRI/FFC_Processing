function [Kout, info] = rs2d_predict_coherent_rf(Kin, phaseOrder, is2D, model)
%RS2D_PREDICT_COHERENT_RF Predict additive RF tones in Cartesian complex k-space.
% K dimensions are [read, centred ky, partition/slice, time, field, receiver].
% is2D treats dimension 3 as independent slices. It defaults to true for
% singleton z and false for multi-partition inputs (legacy 3-D behavior).
% Supported orders include 'auto', 'centric-auto', 'centric-positive',
% 'centric-negative', 'linear', and 'centric2d' (automatic order with a
% 2-D spatial corroboration check for RS2D Centric2D).
% All fits use stored phase-encode ranks. No ADC or frame clock is assumed.
% Processing is deterministic and never modifies Kin. Pass the original data.
% The optional model is a fitted result for this exact source and geometry.
if nargin < 2, phaseOrder = 'auto'; end
phaseOrder = validatestring(phaseOrder, ...
    {'auto','centric-auto','centric-positive','centric-negative','linear','centric2d'});
validateattributes(Kin, {'single','double'}, {'nonempty','nonsparse'});
if nargin<3, is2D = size(Kin,3)==1; end
validateattributes(is2D,{'logical','numeric'},{'scalar','binary'});
is2D = logical(is2D);
shape = size(Kin);
shape(end+1:6) = 1;
nx = shape(1); ny = shape(2); nz = shape(3);
nf = shape(4)*shape(5); nc = shape(6);
Kout = Kin;
info = struct('version',3,'is2D',is2D,'status','no_correction','phaseOrder',phaseOrder, ...
    'shape',shape,'receivers',[]);
if numel(shape)>6 || nx<32 || ny<48 || isreal(Kin)
    info.status = 'unsupported_input';
    return
end
K = reshape(Kin,[nx ny nz nf nc]);
if nargin>=4 && ~isempty(model)
    assert(isequal(model.shape,shape) && model.is2D==is2D && ...
        strcmp(model.phaseOrder,phaseOrder) && model.version==info.version, ...
        'The cached RF model does not match this acquisition layout.');
    info = model;
    Kout = reshape(applyModel(K,info),size(Kin));
    return
end
ky = (0:ny-1).' - floor(ny/2);
distance = abs(ky);
train = distance >= floor(.34*ny);
tune = distance >= floor(.27*ny) & ~train;
audit = distance >= floor(.19*ny) & ~(train | tune);
discoveryFrames = unique(round(linspace(1,nf,min(12,nf))));
discoveryPartitions = unique(round(linspace(1,nz,min(4,nz))));
useSpatialGuard = is2D && strcmp(phaseOrder,'centric2d');
background = [];
if useSpatialGuard
    middle = distance<=floor(.10*ny);
    background = spatialBackground(K,shape(4),middle);
end
for channel = 1:nc
    raw = K(:,:,:,:,channel);
    eligible = all(isfinite(raw(:))) && ~isreal(raw);
    peaks = zeros(0,3); order = phaseOrder; orderScores = [];
    if eligible
        sample = raw(:,train,discoveryPartitions,discoveryFrames);
        [peaks, rank, order, orderScores] = discover(sample, ky, train, phaseOrder);
    end
    rec = struct('phaseOrder',order,'orderScores',orderScores, ...
        'discoveryFrames',discoveryFrames,'discoveryPartitions',discoveryPartitions, ...
        'trainKy',find(train).','tuneKy',find(tune).','auditKy',find(audit).', ...
        'readoutPixel',peaks(:,1)+1,'initialCyclesPerLine',peaks(:,2), ...
        'discoveryScore',peaks(:,3),'tuneReduction',[],'auditReduction',[], ...
        'auditReductionBySide',[],'accepted',[],'reason',{{}},'localSlope',[], ...
        'fittedPartitions',[],'appliedPartitions',[],'coefficient',[],'slope',[], ...
        'readoutCondition',NaN,'centralGuardRows',[], ...
        'centralGuardGain',[],'centralGuardSE',[],'centralGuardKy',[], ...
        'centralSpatialGain',[],'centralSpatialSE',[],'centralSpatialConfirmed',[], ...
        'slice',[],'status','no_correction');
    if ~eligible
        rec.status = 'invalid_receiver';
    elseif ~isempty(peaks)
        requireCenter = is2D || strcmp(order,'linear');
        if is2D
            for z = 1:nz
                planeBackground = [];
                if useSpatialGuard, planeBackground = background(:,:,z); end
                planeRec = predict(raw(:,:,z,:),rec, ...
                    rank, ky, train, tune, audit, peaks, shape(4),requireCenter,planeBackground);
                if any(planeRec.accepted(:))
                    planeRec.status = 'corrected';
                    rec.status = 'corrected';
                end
                if z==1
                    sliceResults = planeRec;
                else
                    sliceResults(z) = planeRec; %#ok<AGROW>
                end
            end
            rec.slice = sliceResults;
        else
            rec = predict(raw,rec,rank,ky,train,tune,audit, ...
                peaks,shape(4),requireCenter,[]);
            if any(rec.accepted(:))
                rec.status = 'corrected';
            end
        end
        if strcmp(rec.status,'corrected'), info.status = 'corrected'; end
    end
    if channel==1
        info.receivers = rec;
    else
        info.receivers(channel) = rec;
    end
end
if all(strcmp({info.receivers.status},'invalid_receiver')), info.status = 'nonfinite_input'; end
Kout = reshape(applyModel(K,info),size(Kin));
end

function [peaks, rank, chosen, scores] = discover(sample, ky, train, order)
orders = {'centric-positive','centric-negative','linear'};
if strcmp(order,'centric-auto')
    orders = orders(1:2);
elseif ~any(strcmp(order,{'auto','centric2d'}))
    orders = {order};
end
nx = size(sample,1); ny = numel(ky);
nRead = 4*2^nextpow2(nx);
nLine = 4*2^nextpow2(ny);
wRead = .5-.5*cos(2*pi*((0:nx-1).'+.5)/nx);
sample = reshape(double(sample),nx,nnz(train),[]);
scores = zeros(1,numel(orders));
bestMerit = -Inf;
for i = 1:numel(orders)
    trialRank = phaseRank(ky,orders{i});
    r = trialRank(train);
    if strcmp(orders{i},'linear')
        % Taper each measured outer band without downweighting it just
        % because a large, unobserved interval lies between the two bands.
        wLine = zeros(size(r));
        negative = ky(train)<0;
        for side = [false true]
            use = negative==side;
            wLine(use) = .5-.5*cos(2*pi*((0:nnz(use)-1)+.5)/nnz(use));
        end
    else
        wLine = sin(pi*(r-min(r)+.5)/(max(r)-min(r)+1)).^2;
    end
    power = zeros(nRead,nLine);
    padded = complex(zeros(nRead,nLine));
    for j = 1:size(sample,3)
        tapered = (sample(:,:,j).*wRead).*wLine.';
        padded(:,r+1) = fftshift(ifft(tapered,nRead,1),1);
        spectrum = abs(fft(padded,[],2)).^2;
        power = power + spectrum/max(sum(spectrum(:)),realmin);
    end
    % Aggregate normalized spectra so anatomy-heavy frames cannot dominate.
    scores(i) = max(power(:))/max(sum(power(:)),realmin);
    if scores(i)>bestMerit
        bestMerit = scores(i);
        rank = trialRank;
        chosen = orders{i};
        bestPower = power;
    end
end
peaks = zeros(0,3);
floorPower = max(median(max(bestPower,[],2)),realmin);
remaining = bestPower;
readGrid = (0:nRead-1).'*nx/nRead;
for t = 1:16
    [power,j] = max(remaining(:));
    score = power/floorPower;
    if score<3.0, break; end
    [x,p] = ind2sub(size(remaining),j);
    qx = mod(x+(-2:0),nRead)+1;
    qp = mod(p+(-2:0),nLine)+1;
    offsetX = peakOffset(log(max(bestPower(qx,p),realmin)));
    offsetP = peakOffset(log(max(bestPower(x,qp),realmin)));
    px = mod((x-1+offsetX)*nx/nRead,nx);
    fp = mod((p-1+offsetP)/nLine,1);
    peaks(end+1,:) = [px,fp,score]; %#ok<AGROW>
    % The model estimates one slow-time carrier per readout tone.
    remaining(abs(mod(readGrid-px+nx/2,nx)-nx/2)<2.25,:) = 0;
end
end

function r = phaseRank(ky, order)
switch order
    case 'centric-positive'
        r = 2*abs(ky) - (ky>0);
    case 'centric-negative'
        r = 2*abs(ky) - (ky<0);
    otherwise
        r = ky-min(ky);
end
% At even matrix sizes the outer negative Nyquist line can be unpaired.
[~, index] = sort(r);
r(index) = 0:numel(ky)-1;
end

function delta = peakOffset(logPower)
a=logPower(1); b=logPower(2); c=logPower(3);
curvature = a-2*b+c;
delta = 0;
if curvature<0, delta = max(-.5,min(.5,.5*(a-c)/curvature)); end
end

function rec = predict(K, rec, r, ky, train, tune, audit, peaks, nt, requireCenter, background)
nx = size(K,1); ny = size(K,2); nz = size(K,3); nf = size(K,4);
np = size(peaks,1);
b = readoutBasis(nx,peaks(:,1)+1);
[U,S,V] = svd(b,'econ');
s = diag(S);
rec.readoutCondition = s(1)/s(end);
if rec.readoutCondition>6
    rec.reason = {'ill_conditioned_readout_basis'};
    return
end
dual = (V./s.')*U';
r = r - mean(r(train));
keep = train | tune | audit;
auditSides = [audit & ky<0, audit & ky>0];
rec.tuneReduction = NaN(np,nf);
rec.auditReduction = NaN(np,nf);
rec.auditReductionBySide = NaN(np,nf,2);
rec.accepted = false(np,nf);
rec.reason = repmat({'weak_or_inconsistent_tuning'},np,nf);
rec.localSlope = false(np,nf);
rec.fittedPartitions = zeros(np,nf);
rec.appliedPartitions = zeros(np,nf);
rec.coefficient = complex(zeros(np,nz,nf));
rec.slope = NaN(np,nz,nf);
bx = fftshift(ifft(ifftshift(b,1),[],1),1);
center = abs(ky)<=floor(.10*ny);
slabs = readoutMargins(K, nt, center);
rec.centralGuardRows = {find(slabs(:,1)),find(slabs(:,2))};
rec.centralGuardKy = find(center).';
rec.centralGuardGain = NaN(np,nf,2);
rec.centralGuardSE = NaN(np,nf,2);
if any(background(:))
    rec.centralSpatialGain = complex(NaN(np,nf));
    rec.centralSpatialSE = NaN(np,nf);
    rec.centralSpatialConfirmed = false(np,nf);
end
canCheckCenter = all(sum(slabs,1)>=6);
rt = r(train);
for f = 1:nf
    block = K(:,:,:,f);
    q = complex(zeros(np,ny,nz));
    selected = double(block(:,keep,:));
    q(:,keep,:) = reshape(dual*reshape(selected,nx,[]),np,nnz(keep),nz);
    predicted = complex(zeros(np,ny,nz));
    for t = 1:np
        y = reshape(q(t,:,:),ny,nz);
        w = peaks(t,2);
        a = mean(exp(-2i*pi*rt*w).*y(train,:),1);
        design = exp(2i*pi*r*w);
        fitted = design*a;
        best = reduction(y,fitted,tune);
        if peaks(t,3)>=20
            % Optimise local slopes on the outer training lines only.
            % Select the resulting family on the tuning band.
            [local, al, wl] = localFit(y(train,:),r,rt,w);
            localScore = reduction(y,local,tune);
            if localScore>best+.0005
                fitted = local;
                w = wl;
                a = al;
                best = localScore;
                rec.localSlope(t,f) = true;
            end
        end
        % Keep only partitions with sufficient coherent evidence in training.
        errorPower = mean(abs(y(train,:)-fitted(train,:)).^2,1);
        coherent = abs(a).^2*nnz(train) >= 8*max(errorPower,realmin);
        fitted(:,~coherent) = 0;
        a(~coherent) = 0;
        rec.fittedPartitions(t,f) = nnz(coherent);
        rec.coefficient(t,:,f) = a;
        rec.slope(t,:,f) = w;
        rec.tuneReduction(t,f) = reduction(y,fitted,tune);
        rec.auditReduction(t,f) = reduction(y,fitted,audit);
        for side = 1:2
            rec.auditReductionBySide(t,f,side) = reduction(y,fitted,auditSides(:,side));
        end
        if rec.tuneReduction(t,f)<.12, continue; end
        rec.reason{t,f} = 'inconsistent_inner_band';
        if rec.auditReduction(t,f)<.08 || any(rec.auditReductionBySide(t,f,:)<-.02)
            continue
        end
        rec.reason{t,f} = 'accepted';
        rec.accepted(t,f) = true;
        predicted(t,:,:) = fitted;
    end
    confirmed = false(np,1);
    if canCheckCenter && any(rec.accepted(:,f))
        [gain, uncertainty] = centralConsistency(block(:,center,:), ...
            predicted(:,center,:), bx, slabs, rec.accepted(:,f));
        rec.centralGuardGain(:,f,:) = reshape(gain,[np 1 2]);
        rec.centralGuardSE(:,f,:) = reshape(uncertainty,[np 1 2]);
        confirmed = any(gain-2.5*uncertainty>.25 & gain+2.5*uncertainty<2,2);
        contradicted = all(gain+3*uncertainty<.20,2);
        rec.reason(contradicted,f) = repmat({'central_extrapolation_contradicted'},nnz(contradicted),1);
        rec.accepted(contradicted,f) = false;
    end
    if any(background(:)) && ~strcmp(rec.phaseOrder,'linear')
        check = rec.accepted(:,f) & ~confirmed;
        if any(check)
            spatialCoefficients = predicted(:,center,:);
            spatialCoefficients(~rec.accepted(:,f),:) = 0;
            [g,se] = centralSpatialConsistency(block(:,center,:),spatialCoefficients, ...
                bx, center, background, check);
            rec.centralSpatialGain(:,f) = g;
            rec.centralSpatialSE(:,f) = se;
            rec.centralSpatialConfirmed(:,f) = abs(g-1)+3*se<.70;
            confirmed = confirmed | rec.centralSpatialConfirmed(:,f);
        end
    end
    if requireCenter
        unconfirmed = rec.accepted(:,f) & ~confirmed;
        rec.reason(unconfirmed,f) = repmat({'central_prediction_unconfirmed'},nnz(unconfirmed),1);
        rec.accepted(unconfirmed,f) = false;
    end
    rec.appliedPartitions(:,f) = rec.fittedPartitions(:,f).*rec.accepted(:,f);
end
end

function out = applyModel(K,info)
out = K;
if ~strcmp(info.status,'corrected'), return; end
for c = 1:numel(info.receivers)
    rec = info.receivers(c);
    if ~strcmp(rec.status,'corrected'), continue; end
    b = readoutBasis(size(K,1),rec.readoutPixel);
    ky = (0:size(K,2)-1).' - floor(size(K,2)/2);
    rank = phaseRank(ky,rec.phaseOrder);
    rank = rank - mean(rank(rec.trainKy));
    if info.is2D
        for z = 1:numel(rec.slice)
            plane = rec.slice(z);
            for f = find(any(plane.accepted,1))
                correction = modelBlock(b,rank,plane,f);
                out(:,:,z,f,c) = K(:,:,z,f,c)-cast(correction,class(K));
            end
        end
    else
        for f = find(any(rec.accepted,1))
            correction = modelBlock(b,rank,rec,f);
            out(:,:,:,f,c) = K(:,:,:,f,c)-cast(correction,class(K));
        end
    end
end
end

function correction = modelBlock(b,r,rec,f)
active = find(rec.accepted(:,f));
nz = size(rec.coefficient,2);
coeff = reshape(rec.coefficient(active,:,f),numel(active),1,nz);
slope = reshape(rec.slope(active,:,f),numel(active),1,nz);
wave = exp(2i*pi*r.' .* slope).*coeff;
correction = reshape(b(:,active)*reshape(wave,numel(active),[]), ...
    size(b,1),numel(r),nz);
end

function b = readoutBasis(nx,pixel)
r = (0:nx-1).' - floor(nx/2);
b = exp(-2i*pi*r*(pixel(:).'-1-floor(nx/2))/nx);
end

function slabs = readoutMargins(K, nt, center)
% A robust late-frame readout profile identifies two disjoint background slabs.
nx = size(K,1);
lastTimes = nt:nt:size(K,4);
indices = unique(round(linspace(1,numel(lastTimes),min(5,numel(lastTimes)))));
profiles = zeros(nx,numel(indices));
for i = 1:numel(indices)
    I = zeros(nx,size(K,2),size(K,3),'like',K);
    I(:,center,:) = K(:,center,:,lastTimes(indices(i)));
    for d = 1:ndims(I), I = ifft(ifftshift(I,d),[],d); end
    profile = sort(reshape(abs(fftshift(I,1)),nx,[]),2);
    profiles(:,i) = profile(:,ceil(.80*size(profile,2)));
end
raw = median(profiles,2);
profile = median([circshift(raw,-3) circshift(raw,-2) circshift(raw,-1) ...
    raw circshift(raw,1) circshift(raw,2) circshift(raw,3)],2);
ordered = sort(profile);
reference = ordered(ceil(.90*nx));
x = (0:nx-1).'/nx;
trusted = profile < .20*reference;
slabs = [trusted & x<.23, trusted & x>=.77];
end

function [gain, uncertainty] = centralConsistency(K, coeff, basis, slabs, active)
% Test inner ky with the readout-image samples in each background slab.
% Use ky-clustered residuals to allow correlation across partitions and x.
nx = size(K,1); nky = size(K,2); nz = size(K,3);
nt = size(coeff,1);
C = reshape(coeff,nt,[]);
Y = fftshift(ifft(ifftshift(double(reshape(K,nx,[])),1),[],1),1);
Y = Y - basis*C;
gain = NaN(nt,2); uncertainty = NaN(nt,2);
for side = 1:2
    xs = slabs(:,side);
    for t = find(active).'
        estimate = basis(xs,t)*C(t,:);
        perKy = reshape(sum(conj(estimate).*Y(xs,:),1),nky,nz);
        estEnergy = reshape(sum(abs(estimate).^2,1),nky,nz);
        weight = sum(estEnergy,2);
        total = sum(weight);
        if total<=realmin, continue; end
        residual = real(sum(perKy,2));
        delta = sum(residual)/total;
        gain(t,side) = 1+delta;
        uncertainty(t,side) = sqrt(nky/(nky-1)*sum((residual-delta*weight).^2))/total;
    end
end
end

function background = spatialBackground(K, nt, center)
% Coil-combined central-only support, with a robust reference across fields.
nx = size(K,1); ny = size(K,2); nz = size(K,3);
lastTimes = nt:nt:size(K,4);
ids = unique(round(linspace(1,numel(lastTimes),min(5,numel(lastTimes)))));
volumes = zeros(nx,ny,nz,numel(ids),'single');
for i = 1:numel(ids)
    central = zeros(nx,ny,nz,1,size(K,5),'like',K);
    central(:,center,:,:,:) = K(:,center,:,lastTimes(ids(i)),:);
    central(~isfinite(central)) = 0;
    image = centeredInverse2(central);
    volumes(:,:,:,i) = sqrt(sum(abs(image).^2,5));
end
ref = median(volumes,4);
% A readout median suppresses narrow RF stripes in the support reference.
% All regions around uncertain tissue are excluded from the test.
ref = median(cat(4,ref,circshift(ref,1,1),circshift(ref,-1,1), ...
    circshift(ref,2,1),circshift(ref,-2,1)),4);
background = false(nx,ny,nz);
for z = 1:nz
    p = ref(:,:,z);
    smoothed = zeros(nx,ny);
    for dx = -2:2
        for dy = -2:2
            smoothed = smoothed + double(circshift(p,[dx dy]));
        end
    end
    levels = sort(smoothed(:));
    upper = levels(ceil(.95*numel(levels)));
    noiseFloor = levels(ceil(.20*numel(levels)));
    if upper<=0 || noiseFloor>=.4*upper, continue; end
    threshold = max(.13*upper, 2.5*noiseFloor);
    support = smoothed>=threshold;
    excluded = support;
    for dx = [-4 0 4]
        for dy = [-6 0 6]
            excluded = excluded | circshift(support,[dx dy]);
        end
    end
    background(:,:,z) = ~excluded;
end
end

function [gain, se] = centralSpatialConsistency(K, coeff, basis, center, background, active)
% Test the predicted complex template on off-tissue central-line images.
% The short ky bandwidth correlates image samples in y. The SE includes
% a bandwidth penalty so the zero-filled y samples cannot masquerade as
% independent measurements.
nx = size(K,1); ny = numel(center);
hybrid = complex(zeros(nx,ny));
hybrid(:,center) = fftshift(ifft(ifftshift(double(K),1),[],1),1);
target = fftshift(ifft(ifftshift(hybrid,2),[],2),2);
wave = complex(zeros(size(coeff,1),ny));
wave(:,center) = coeff;
wave = fftshift(ifft(ifftshift(wave,2),[],2),2);
residual = target-basis*wave;
gain = complex(NaN(size(active)));
se = NaN(size(active));
residual = residual(background);
sigma = sqrt(median(abs(residual).^2)/log(2));
for t = find(active).'
    estimate = basis(:,t)*wave(t,:);
    on = sum(abs(estimate(background)).^2);
    total = sum(abs(estimate(:)).^2);
    if on < .12*total || on<=realmin, continue; end
    gain(t) = 1+sum(conj(estimate(background)).*residual)/on;
    se(t) = sigma*sqrt(ny/nnz(center)/on);
end
end

function out = centeredInverse2(in)
out = fftshift(ifft(ifftshift(in,1),[],1),1);
out = fftshift(ifft(ifftshift(out,2),[],2),2);
end

function [pred, a, w] = localFit(y, r, rt, seed)
nz = size(y,2);
w = zeros(1,nz);
grid = seed + (-8:8)*.0005;
coarse = abs(exp(-2i*pi*grid(:)*rt.')*y).^2;
[~,idx] = max(coarse,[],1);
for z = 1:nz
    j = idx(z);
    offset = 0;
    if j>1 && j<numel(grid)
        offset = peakOffset(log(max(coarse(j-1:j+1,z),realmin)));
    end
    w(z) = grid(j)+offset*.0005;
end
a = mean(exp(-2i*pi*rt*w).*y,1);
pred = exp(2i*pi*r*w).*a;
end

function value = reduction(y, yhat, rows)
original = y(rows,:);
residual = original-yhat(rows,:);
value = (sum(abs(original(:)).^2)-sum(abs(residual(:)).^2))/ ...
    max(sum(abs(original(:)).^2),realmin);
end
