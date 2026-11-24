function [kspaceCorrected, info] = rs2d_remove_coherent_tones(kspace, varargin)
%RS2D_REMOVE_COHERENT_TONES v7.1: adaptive residual-tone refitting.
%
%   [Kcorr,info] = rs2d_remove_coherent_tones(K, ...
%       'SpatialDimensions',3,'CoilDimension',0, ...
%       'PhaseOrder','centric-alternating','KZFast',false, ...
%       'HarmonicMode','auto');
%
% Runs the proven v5 temporal-tone estimator, then selects UNDER-SUBTRACTED
% tones and independently calibrates their remaining complex amplitude from
% high-confidence OUTSIDE-ANATOMY ghost samples. The SAME fitted additive
% tone model is subtracted throughout 3D space, including inside tissue.
% This differs from v6, which suppressed subtraction inside the tissue mask.
%
% Existing v5 name-value arguments are forwarded unchanged to v5core.
% New options (name/value):
%  HarmonicMode                 'auto' (3D), 'off'
%  HarmonicReadoutPositions     [] 1-based readout pixels; e.g. 64.516314
%  MaxHarmonicCandidates        3, max candidate models per frame
%  MaxHarmonicsAccepted         1, max extra subtractions per frame
%  HarmonicBackgroundROI        [] logical [nx ny nz] verified ghost-only ROI
%  HarmonicSupportMask          [] verified brain support (recommended)
%  HarmonicModelKspace          [] EXTERNAL fitted additive residual model
%  BaselineKspace               [] already corrected v5 data, same shape
%  BaselineInfo                 [] v5 information (if baseline supplied)
%  MinHarmonicHoldoutReduction  0.70, held-out complex background reduction
%  MaxHarmonicGain              1.35; fit |alpha| must be <= this
%  MaxHarmonicPhase             0.30 rad
%  MaxHarmonicAnatomyChange     0.05, versus ORIGINAL K, whole tissue
%  MaxHarmonicPerSliceChange    0.10, versus ORIGINAL K, each slice
%  MinHarmonicStrength          0.25; lower accepted scale -> skip
%  HarmonicAutoPercentile       99, choose predicted-tone hotspots in bg
%  HarmonicResidualRefit        false (experimental); refit tones from v5 residual
%  HarmonicVerbose              false
%
% Input/output same shape and class. This is RESEARCH/EXPERIMENTAL and its
% background holdout is spatially correlated, NOT proof of tissue recovery.
% Inspect original/corrected T1 and R1 curves before quantitative use.
% No Image Processing Toolbox required. No external files are read.
%
% 'HarmonicModelKspace' + 'BaselineKspace' reproduce the manual target
% experiment WITHOUT re-fitting; automatically fitted mode uses v5 toneTests.

hNames={'harmonicmode','harmonicreadoutpositions','maxharmoniccandidates', ...
 'maxharmonicsaccepted','harmonicbackgroundroi','harmonicsupportmask', ...
 'harmonicmodelkspace','baselinekspace','baselineinfo', ...
 'minharmonicholdoutreduction','maxharmonicgain','maxharmonicphase', ...
 'maxharmonicanatomychange','maxharmonicperslicechange', ...
 'minharmonicstrength','harmonicautopercentile','harmonicverbose', ...
 'harmonicresidualrefit'};
if mod(numel(varargin),2)~=0
    error('rs2d_remove_coherent_tones:NameValue','Expected name/value pairs.');
end
legacy={};harmonic={};
for ii=1:2:numel(varargin)
    key=varargin{ii};
    if ~(ischar(key)||isstring(key))
        error('rs2d_remove_coherent_tones:NameValue','Option names must be text.');
    end
    if ismember(lower(char(key)),hNames)
        harmonic(end+1:end+2)=varargin(ii:ii+1); %#ok<AGROW>
    else
        legacy(end+1:end+2)=varargin(ii:ii+1); %#ok<AGROW>
    end
end
p=inputParser;
addParameter(p,'HarmonicMode','auto',@(x)ischar(x)||isstring(x));
addParameter(p,'HarmonicReadoutPositions',[],@(x)isnumeric(x)&&all(isfinite(x(:))));
addParameter(p,'MaxHarmonicCandidates',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==fix(x));
addParameter(p,'MaxHarmonicsAccepted',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x==fix(x));
addParameter(p,'HarmonicBackgroundROI',[],@(x)isempty(x)||islogical(x)||isnumeric(x));
addParameter(p,'HarmonicSupportMask',[],@(x)isempty(x)||islogical(x)||isnumeric(x));
addParameter(p,'HarmonicModelKspace',[],@(x)isempty(x)||isnumeric(x));
addParameter(p,'BaselineKspace',[],@(x)isempty(x)||isnumeric(x));
addParameter(p,'BaselineInfo',[],@(x)isempty(x)||isstruct(x));
addParameter(p,'MinHarmonicHoldoutReduction',0.70,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MaxHarmonicGain',1.35,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MaxHarmonicPhase',0.30,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=pi);
addParameter(p,'MaxHarmonicAnatomyChange',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MaxHarmonicPerSliceChange',0.10,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MinHarmonicStrength',0.25,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'HarmonicAutoPercentile',99,@(x)isnumeric(x)&&isscalar(x)&&x>=90&&x<=99.9);
addParameter(p,'HarmonicVerbose',false,@(x)islogical(x)||isscalar(x));
addParameter(p,'HarmonicResidualRefit',false,@(x)islogical(x)||isscalar(x));
parse(p,harmonic{:});h=p.Results;
mode=lower(strtrim(char(h.HarmonicMode)));
if ~ismember(mode,{'auto','off'})
    error('rs2d_remove_coherent_tones:HarmonicMode','Use auto or off.');
end

if isempty(h.BaselineKspace)
    [Kbase,baseInfo]=v5core(kspace,legacy{:});
else
    if ~isequal(size(h.BaselineKspace),size(kspace))
        error('rs2d_remove_coherent_tones:BaselineSize','BaselineKspace must match kspace.');
    end
    Kbase=cast(h.BaselineKspace,'like',kspace);
    if ~isempty(h.BaselineInfo)
        baseInfo=h.BaselineInfo;
    else
        baseInfo=struct('status','supplied','toneTests',[]);
    end
end
kspaceCorrected=Kbase;
info=baseInfo;info.harmonicVersion='7.1';
info.harmonicResidualRefit=logical(h.HarmonicResidualRefit);
info.harmonicStatus='off';info.harmonicTests=[];info.harmonicNumAccepted=0;
info.harmonicFrameSummary=[];info.harmonicSupportMask=[];
if strcmp(mode,'off')
    return;
end
D=legacyOption(legacy,'SpatialDimensions',2);
if D~=3
    info.harmonicStatus='skipped_non3d';return;
end
if ~isfield(baseInfo,'status') || ...
    (~ismember(baseInfo.status,{'corrected','no_correction','supplied'}))
    info.harmonicStatus='skipped_baseline_status';return;
end
originalShape=size(kspace);nd=ndims(kspace);
nx=size(kspace,1);ny=size(kspace,2);nz=size(kspace,3);
if nx<16||ny<16||nz<2
    info.harmonicStatus='skipped_small_input';return;
end
cdim=legacyOption(legacy,'CoilDimension',[]);
if isempty(cdim)
    if nd>3, cdim=nd; else, cdim=0; end
end
if cdim>nd || (cdim>0 && cdim<=3)
    error('rs2d_remove_coherent_tones:CoilDimension','Invalid coil dimension.');
end
if cdim==0
    perm=1:nd;extra=4:nd;nc=1;
else
    extra=setdiff(4:nd,cdim,'stable');perm=[1:3 extra cdim];nc=size(kspace,cdim);
end
if isempty(extra),nf=1;else,nf=prod(originalShape(extra));end
Korig=reshape(double(permute(kspace,perm)),[nx ny nz nf nc]);
KbaseR=reshape(double(permute(Kbase,perm)),[nx ny nz nf nc]);
if ~isempty(h.HarmonicModelKspace)
    if ~isequal(size(h.HarmonicModelKspace),size(kspace))
        error('rs2d_remove_coherent_tones:ModelSize','HarmonicModelKspace must match input.');
    end
    models=reshape(double(permute(h.HarmonicModelKspace,perm)),[nx ny nz nf nc]);
else
    models=[];
end
if ~isempty(h.HarmonicSupportMask)
    if numel(h.HarmonicSupportMask)~=nx*ny*nz
        error('rs2d_remove_coherent_tones:SupportSize','HarmonicSupportMask must be [kx ky kz].');
    end
    support=reshape(logical(h.HarmonicSupportMask),[nx ny nz]);
else
    sf=unique(max(1,min(nf,round(linspace(1,nf,min(5,nf))))));
    support=makeHarmonicSupport(Korig,sf);
    support=largest3DComponent(support);
end
if mean(support(:))>.80 || nnz(support)<150
    info.harmonicStatus='skipped_unreliable_support';return;
end
info.harmonicSupportMask=support;
if ~isempty(h.HarmonicBackgroundROI)
    if numel(h.HarmonicBackgroundROI)~=nx*ny*nz
        error('rs2d_remove_coherent_tones:ROI','HarmonicBackgroundROI must be [kx ky kz].');
    end
    userROI=reshape(logical(h.HarmonicBackgroundROI),[nx ny nz]);
    if any(userROI(:)&support(:))
        error('rs2d_remove_coherent_tones:ROIOverlap','Background ROI overlaps anatomy mask.');
    end
else
    userROI=[];
end
% v5 already checked complex validity; don't alter on missing candidate.
work=KbaseR;
recordTemplate=struct('frame',0,'readoutPosition',NaN,'status','',...
 'complexGain',complex(NaN,NaN),'holdoutReduction',NaN, ...
 'backgroundVoxelCount',0,'strength',0,'anatomyChange',NaN, ...
 'maxSliceChange',NaN,'backgroundReduction',NaN);
tests=repmat(recordTemplate,max(1,nf*h.MaxHarmonicCandidates),1);nTests=0;
summary=zeros(nf,5);numAccepted=0;
for f=1:nf
    kc=reshape(work(:,:,:,f,:),[nx ny nz nc]);
    kref=reshape(Korig(:,:,:,f,:),[nx ny nz nc]);
    imgRef=complexImages3D(kref);img=complexImages3D(kc);
    tissueRef=sqrt(sum(abs(imgRef).^2,4));
    tissue0=sqrt(sum(abs(img).^2,4));
    denom=max(mean(tissueRef(support)),eps);
    sliceDen=ones(nz,1);
    for iz=1:nz
        zz=support(:,:,iz);pp=tissueRef(:,:,iz);
        if any(zz(:)),sliceDen(iz)=max(mean(pp(zz)),eps);end
    end
    originalBackground=sum(tissue0(~support).^2);
    candidateModels={}; candidatePositions=[]; candidateScores=[];
    if ~isempty(models)
        candidateModels={reshape(models(:,:,:,f,:),[nx ny nz nc])};
        candidatePositions=NaN;candidateScores=Inf;
    else
        % Reconstruct residual complex tone directly from v5 saved fit.
        % This avoids re-estimating it with any anatomical pixels.
        if isfield(baseInfo,'toneTests') && ~isempty(baseInfo.toneTests) && ...
                isfield(baseInfo,'acquisitionIndex') && ...
                isfield(baseInfo,'envelopeKnotTimes')
            for t=1:numel(baseInfo.toneTests)
                tt=baseInfo.toneTests(t);
                if tt.frame~=f || isempty(tt.envelopeCoefficients) || ...
                   ~isfinite(tt.phaseCyclesPerReadout) || tt.strength>=0.98 || ...
                   tt.heldOutImprovement<0.20 || tt.trainingCoherence<0.06
                    continue;
                end
                if ~isempty(h.HarmonicReadoutPositions) && ...
                   all(abs(tt.readoutPosition-h.HarmonicReadoutPositions)>1.1)
                   continue;
                end
                estimated=residualModelFromTest(tt,baseInfo,nx,ny,nz,nc);
                if isempty(estimated),continue;end
                candidateModels{end+1}=estimated; %#ok<AGROW>
                candidatePositions(end+1)=tt.readoutPosition; %#ok<AGROW>
                candidateScores(end+1)=tt.heldOutImprovement* ...
                    (1-tt.strength)*max(tt.trainingCoherence,0); %#ok<AGROW>
            end
        end
        % If a v5 baseline was supplied without its fit diagnostics, allow
        % explicit positions to be refitted on the residual, unrestricted.
        if isempty(candidateModels) && ~isempty(h.HarmonicReadoutPositions)
            for j=1:min(h.MaxHarmonicCandidates,numel(h.HarmonicReadoutPositions))
                x=h.HarmonicReadoutPositions(j);
                if ~isfinite(x)||x<1||x>nx,continue;end
                [estimated,ok]=fitResidualV5(kc,x,support,legacy);
                if ok
                    candidateModels{end+1}=estimated; %#ok<AGROW>
                    candidatePositions(end+1)=x; %#ok<AGROW>
                    candidateScores(end+1)=1; %#ok<AGROW>
                end
            end
        end
    end
    % Optional v7.1 fallback: refit the REMAINING temporal tones directly
    % from Kbase. Unlike v7, this searches even if the original v5 tone
    % was fully accepted or not among its partly clipped candidates.
    % Each model still needs independent background holdout and tissue QC.
    if h.HarmonicResidualRefit
        if ~isempty(h.HarmonicReadoutPositions)
            refitPositions=double(h.HarmonicReadoutPositions(:)).';
        elseif isfield(baseInfo,'refinedReadoutPositions')
            refitPositions=double(baseInfo.refinedReadoutPositions(:)).';
        else
            refitPositions=[];
        end
        for jj=1:numel(refitPositions)
            x=refitPositions(jj);
            if ~isfinite(x)||x<1||x>nx,continue;end
            [estimated,ok]=fitResidualV5(kc,x,support,legacy);
            if ~ok,continue;end
            candidateModels{end+1}=estimated; %#ok<AGROW>
            candidatePositions(end+1)=x; %#ok<AGROW>
            % Prefer direct residual fits over older, scaled v5 models.
            candidateScores(end+1)=2-1e-5*jj; %#ok<AGROW>
        end
    end
    if ~isempty(candidateScores)
        [~,ix]=sort(candidateScores,'descend');
        ix=ix(1:min(numel(ix),h.MaxHarmonicCandidates));
        candidateModels=candidateModels(ix);candidatePositions=candidatePositions(ix);
    end
    naccepted=0;
    for j=1:numel(candidateModels)
        if naccepted>=h.MaxHarmonicsAccepted,break;end
        modk=candidateModels{j};
        modimg=complexImages3D(modk);
        % Fit one shared complex alpha across receiver channels using only
        % background ghost; anatomical support NEVER enters the alpha fit.
        if isempty(userROI)
            roi=automaticBackgroundROI(img,modimg,support,h.HarmonicAutoPercentile);
        else
            roi=userROI;
        end
        nTests=nTests+1;rec=recordTemplate;
        rec.frame=f;rec.readoutPosition=candidatePositions(j);
        rec.backgroundVoxelCount=nnz(roi);
        if nnz(roi)<64
            rec.status='too_few_background_pixels';tests(nTests)=rec;continue;
        end
        [a,b,c]=ndgrid(1:nx,1:ny,1:nz);
        train=roi & mod(b+c,2)==0;
        valid=roi & ~train;
        if nnz(train)<24||nnz(valid)<24
            rec.status='too_few_holdout_pixels';tests(nTests)=rec;continue;
        end
        % Index all coils without a loop in the fit, but keep the same
        % scalar complex alpha for every channel of the frame.
        in=reshape(img,[],nc); no=reshape(modimg,[],nc);
        tr=find(train(:));va=find(valid(:));
        aTrain=no(tr,:);bTrain=in(tr,:);
        den=sum(abs(aTrain(:)).^2);
        if den<eps
            rec.status='no_model_energy';tests(nTests)=rec;continue;
        end
        alpha=sum(conj(aTrain(:)).*bTrain(:))/den;
        rec.complexGain=alpha;
        if abs(alpha)>h.MaxHarmonicGain || abs(angle(alpha))>h.MaxHarmonicPhase
            rec.status='gain_out_of_bounds';tests(nTests)=rec;continue;
        end
        aVal=no(va,:);bVal=in(va,:);
        energy0=sum(abs(bVal(:)).^2);
        frac=1-sum(abs(bVal(:)-alpha*aVal(:)).^2)/max(energy0,eps);
        rec.holdoutReduction=frac;
        if frac<h.MinHarmonicHoldoutReduction
            rec.status='failed_background_holdout';tests(nTests)=rec;continue;
        end
        % Enforce tissue change w.r.t. ORIGINAL K before permitting full
        % subtraction. Scale down only if necessary, rather than weakening
        % the interference estimate in background before calibration.
        [ok,change,local]=harmonicLimits(img-alpha*modimg,tissueRef,support, ...
            denom,sliceDen,h.MaxHarmonicAnatomyChange,h.MaxHarmonicPerSliceChange);
        if ok
            strength=1;
        else
            lo=0;hi=1;
            for it=1:12
                md=(lo+hi)*.5;
                [safe,~,~]=harmonicLimits(img-md*alpha*modimg,tissueRef,support, ...
                    denom,sliceDen,h.MaxHarmonicAnatomyChange,h.MaxHarmonicPerSliceChange);
                if safe,lo=md;else,hi=md;end
            end
            strength=lo;
        end
        if strength<h.MinHarmonicStrength
            rec.status='anatomy_change_limit';tests(nTests)=rec;continue;
        end
        img=img-strength*alpha*modimg;
        kc=fft3DCoils(img);
        numAccepted=numAccepted+1;naccepted=naccepted+1;
        [~,change,local]=harmonicLimits(img,tissueRef,support, ...
            denom,sliceDen,Inf,Inf);
        mnew=sqrt(sum(abs(img).^2,4));
        rec.status='corrected';rec.strength=strength;
        rec.anatomyChange=change;rec.maxSliceChange=max(local);
        rec.backgroundReduction=1-sum(mnew(~support).^2)/max(originalBackground,eps);
        tests(nTests)=rec;
        if h.HarmonicVerbose
            fprintf('v7 f%d tone %.3f: alpha=%+.4f%+.4fi, holdout %.1f%%, strength %.3f, anatomy %.2f%%\n',...
                f,rec.readoutPosition,real(alpha),imag(alpha),100*frac,strength,100*change);
        end
    end
    work(:,:,:,f,:)=reshape(kc,[nx ny nz 1 nc]);
    newMag=sqrt(sum(abs(img).^2,4));
    [~,e,el]=harmonicLimits(img,tissueRef,support,denom,sliceDen,Inf,Inf);
    summary(f,:)=[f,naccepted,1-sum(newMag(~support).^2)/max(originalBackground,eps),e,max(el)];
end
if cdim==0
    restored=reshape(work,originalShape);
else
    restored=ipermute(reshape(work,originalShape(perm)),perm);
end
kspaceCorrected=cast(restored,'like',kspace);
info.harmonicTests=tests(1:nTests);
info.harmonicNumAccepted=numAccepted;
info.harmonicFrameSummary=summary;
if numAccepted>0
    info.harmonicStatus='corrected';
else
    info.harmonicStatus='no_correction';
end
end

function val=legacyOption(options,name,default)
val=default;
for j=1:2:numel(options)
    if strcmpi(char(options{j}),name),val=options{j+1};end
end
end

function model=residualModelFromTest(tt,baseInfo,nx,ny,nz,nc)
model=[];
coef=double(tt.envelopeCoefficients);
knots=double(baseInfo.envelopeKnotTimes(:).');
if isempty(coef)||numel(knots)~=size(coef,1)||size(coef,2)~=nc,return;end
pidx=double(baseInfo.acquisitionIndex(:));
if numel(pidx)~=ny*nz,return;end
if size(coef,1)==1
    env=repmat(coef,ny*nz,1);
else
    B=interp1(knots,eye(numel(knots)),pidx,'linear');
    env=B*coef;
end
phase=exp(2i*pi*tt.phaseCyclesPerReadout*pidx);
pred=bsxfun(@times,phase,env);
bn=(0:nx-1)'-floor(nx/2);
x=double(tt.readoutPosition);
bx=exp(-2i*pi*(x-(floor(nx/2)+1))*bn/nx);
model=bsxfun(@times,reshape(bx,[nx 1 1 1]), ...
    reshape(pred,[1 ny nz nc]));
model=(1-double(tt.strength))*model;
end

function [model,ok]=fitResidualV5(k,x,support,legacy)
model=[];ok=false;
[nx,ny,nz,nc]=size(k);
try
    args={'SpatialDimensions',3,'CoilDimension',iff(nc==1,0,4), ...
      'ToneReadoutPositions',x,'RefineReadout',false,'MaxCandidates',1,...
      'MaxTones',1,'EnvelopeKnots',legacyOption(legacy,'EnvelopeKnots',24),...
      'EnvelopeRidge',legacyOption(legacy,'EnvelopeRidge',5), ...
      'MaxTotalAnatomyChange',1,'MaxPerSliceAnatomyChange',1, ...
      'SupportMask',support,'MinValidationReduction',0.20};
    if nc==1
        kInput=reshape(k,[nx ny nz]);
    else
        kInput=k;
    end
    [kout,i]=v5core(kInput,args{:});
    if strcmp(i.status,'corrected')
        model=k-reshape(double(kout),[nx ny nz nc]);
        ok=true;
    end
catch
    ok=false;
end
end

function v=iff(tf,a,b)
if tf,v=a;else,v=b;end
end

function roi=automaticBackgroundROI(im,modimg,support,pct)
a=sqrt(sum(abs(im).^2,4));
n=sqrt(sum(abs(modimg).^2,4));
background=~support;
roi=false(size(support));
if nnz(background)<100,return;end
baseMed=median(a(background));
for q=[pct,97.5,95]
    th=localQuantile(n(background),q/100);
    test=background & n>=th & a>max(2.5*baseMed,eps);
    if nnz(test)>=64
        roi=test;return;
    end
end
end

function [ok,globalChange,local]=harmonicLimits(img,ref,mask,den,sliceDen,maxGlobal,maxSlice)
m=sqrt(sum(abs(img).^2,4));
globalChange=mean(abs(m(mask)-ref(mask)))/max(den,eps);
nz=size(m,3);local=zeros(1,nz);
for z=1:nz
    mm=m(:,:,z); rr=ref(:,:,z); zz=mask(:,:,z);
    if any(zz(:))
        local(z)=mean(abs(mm(zz)-rr(zz)))/max(sliceDen(z),eps);
    end
end
ok=globalChange<=maxGlobal+1e-9 && max(local)<=maxSlice+1e-9;
end

function m=complexImages3D(k)
for d=1:3,k=ifftshift(k,d);end
for d=1:3,k=ifft(k,[],d);end
for d=1:3,k=fftshift(k,d);end
m=k;
end
function k=fft3DCoils(im)
k=im;
for d=1:3,k=ifftshift(k,d);end
for d=1:3,k=fft(k,[],d);end
for d=1:3,k=fftshift(k,d);end
end

function support=makeHarmonicSupport(K,sf)
% Mask for background-calibrated subtraction must NOT include persistent
% bright ghosts: deliberately tighter threshold than v5's loose support.
[nx,ny,nz,~,nc]=size(K);
a=zeros(nx,ny,nz);
for f=sf
    k=reshape(K(:,:,:,f,:),[nx ny nz nc]);
    m=imageMagnitude(k,3);
    a=a+m/max(localQuantile(m(:),.995),eps)/numel(sf);
end
a=movmedian(a,11,1);
a=movmedian(a,5,2);
a=gaussianSmooth(a,1,2.5);
a=gaussianSmooth(a,2,2.0);
a=gaussianSmooth(a,3,0.6);
% The lower (.18) v5 threshold dilates the support into the lower ghost.
% .26 is chosen from the existing 3D brain reference dataset, NOT universal.
support=a>.26*localQuantile(a(:),.995);
end

function out=largest3DComponent(mask)
% Simple base-MATLAB 6-connected component; keeps only main anatomy.
[nx,ny,nz]=size(mask);
ids=find(mask);vis=false(size(mask));
q=zeros(numel(ids),1);best=[];
for ii=1:numel(ids)
    root=ids(ii);if vis(root),continue;end
    head=1;tail=1;q(1)=root;vis(root)=true;
    while head<=tail
        id=q(head);head=head+1;
        [x,y,z]=ind2sub([nx ny nz],id);
        nei=[x-1 y z; x+1 y z; x y-1 z; x y+1 z; x y z-1; x y z+1];
        for jj=1:6
            xp=nei(jj,1);yp=nei(jj,2);zp=nei(jj,3);
            if xp<1||xp>nx||yp<1||yp>ny||zp<1||zp>nz,continue;end
            loc=sub2ind([nx ny nz],xp,yp,zp);
            if mask(loc)&&~vis(loc)
                tail=tail+1;q(tail)=loc;vis(loc)=true;
            end
        end
    end
    if tail>numel(best),best=q(1:tail);end
end
out=false(size(mask));out(best)=true;
if nnz(out)<.05*numel(mask),out=mask;end
end

% ================== Unmodified v5 temporal-tone engine ==================

function [kspaceCorrected, info] = v5core(kspace, varargin)
%RS2D_REMOVE_COHERENT_TONES v5.0 -- slice-protected local-envelope tone subtraction.
%
%   [Kcorr, info] = rs2d_remove_coherent_tones(K, ...
%      'SpatialDimensions', 3, 'CoilDimension', 0, ...
%      'PhaseOrder', 'centric-alternating', 'KZFast', false);
%
% INPUT: COMPLEX, CENTRED, Cartesian k-space [kx ky (kz) otherDims ...].
% READOUT = dim 1; KY = dim 2; KZ = dim 3 for SpatialDimensions=3.
% Remaining dimensions are independent FCI frames, except the explicit coil
% dimension. For 5D FCI [kx ky kz evolutionTime evolutionField], ALWAYS set
% 'SpatialDimensions',3,'CoilDimension',0.
%
% OUTPUT: same shape and numeric class as input; second output contains
% fitted tone properties, train/hold-out improvements and image-change metrics.
% No external files, toolboxes, or MATLAB global state required.
%
% Compared with v4: (1) more local 24-knot complex envelopes instead of
% 8 knots, (2) independent safeguards for every reconstructed 3D slice.
% Shared tone candidates, train-only fitting and held-out validation remain.
% v5 additionally constrains mean anatomy change within EACH reconstructed
% kz slice. It does not model rapid frequency drift unless separately fitted.
% Temporal frequency/phase increment can differ between frames. Calibration
% frequency priors can be used when relevant RF/ADC settings are known.
%
% CAVEAT: Fitting in the presence of MR anatomy is an inverse problem.
% Retain the uncorrected data and verify relaxometry measurements separately.
% No model can restore information irreversibly lost through signal averaging.
%
% Options (defaults):
%   SpatialDimensions          2 (explicit 3 for 3D)
%   CoilDimension             [] auto last axis; 0 none; otherwise index
%   PhaseOrder                auto: linear 2D, centric-alternating 3D
%   KZFast                    false (ky complete before next kz)
%   AcquisitionIndex          [] or [ny nz] integer permutation 0:P-1
%   ToneReadoutPositions      [] automatic; else MATLAB 1-based x positions
%   CandidateFrames           5 representative frames for shared dictionary
%   MaxCandidates             16 (tested even if earlier candidates fail)
%   MaxTones                  8 (maximum ACCEPTED corrections per frame)
%   MinCandidateCoherence     0.14
%   MaxThetaSpread            0.016 cycles/readout
%   RefineReadout             true; sub-bin x-frequency refinement
%   EnvelopeKnots             24; piecewise-linear knots over P readouts
%   EnvelopeRidge             5; penalise envelope second differences
%   MinTemporalCoherence      0.06
%   MinValidationReduction    0.20 (strictly held-out readouts)
%   ValidationFold            3 of 0..3 (25%% held-out)
%   MaxTotalAnatomyChange     0.02 relative mean magnitude difference
%   MaxPerSliceAnatomyChange  0.06; maximum mean magnitude change for any slice
%   MinCorrectionStrength     0.15
%   SupportMask               [] auto, or logical [nx ny (nz)]
%   ReferenceFrequenciesHz    [] or known FID tone frequencies, optional
%   ReferenceReceiverHz       required with ReferenceFrequenciesHz
%   ImagingReceiverHz         required with ReferenceFrequenciesHz
%   DwellSeconds              required with ReferenceFrequenciesHz
%   FrequencySign             +/-1, check demodulation convention
%   Verbose                   false
%   Plot                      false
%   ErrorOnSmallInput         false
%
% If ReferenceFrequenciesHz are given, they are converted to additional
% readout-position candidates. They DO NOT supply interference amplitude
% or phase, nor do they prove the source remains stable during imaging.
%
% Reference implementation tested in Python on the uploaded brain dataset.
% MATLAB port and synthetic regression test supplied; no MATLAB runtime
% was available in the authoring environment for execution.

validateattributes(kspace, {'single','double'}, {'nonempty'}, mfilename, 'kspace', 1);
if isreal(kspace)
    error('rs2d_remove_coherent_tones:RealInput', 'Expected complex k-space.');
end
p = inputParser;
p.FunctionName = mfilename; p.KeepUnmatched = true;
addParameter(p,'SpatialDimensions',2,@(x) isequal(x,2)||isequal(x,3));
addParameter(p,'CoilDimension',[],@(x) isempty(x)||(isscalar(x)&&isnumeric(x)&&x>=0&&x==fix(x)));
addParameter(p,'PhaseOrder','auto',@(x) ischar(x)||isstring(x));
addParameter(p,'KZFast',false,@(x) islogical(x)||isscalar(x));
addParameter(p,'AcquisitionIndex',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'ToneReadoutPositions',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'CandidateFrames',5,@(x) isscalar(x)&&x>=1&&x==fix(x));
addParameter(p,'MaxCandidates',16,@(x) isscalar(x)&&x>=1&&x==fix(x));
addParameter(p,'MaxTones',8,@(x) isscalar(x)&&x>=0&&x==fix(x));
addParameter(p,'MinCandidateCoherence',0.14,@(x) isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MaxThetaSpread',0.016,@(x) isscalar(x)&&x>0&&x<=0.5);
addParameter(p,'RefineReadout',true,@(x) islogical(x)||isscalar(x));
addParameter(p,'EnvelopeKnots',24,@(x) isscalar(x)&&x>=1&&x==fix(x));
addParameter(p,'EnvelopeRidge',5,@(x) isscalar(x)&&x>=0);
addParameter(p,'MinTemporalCoherence',0.06,@(x) isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MinValidationReduction',0.20,@(x) isscalar(x)&&x>=0&&x<=1);
addParameter(p,'ValidationFold',3,@(x) isscalar(x)&&ismember(x,0:3));
addParameter(p,'MaxTotalAnatomyChange',0.02,@(x) isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MaxPerSliceAnatomyChange',0.06,@(x) isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MinCorrectionStrength',0.15,@(x) isscalar(x)&&x>=0&&x<=1);
addParameter(p,'SupportMask',[],@(x) isempty(x)||isnumeric(x)||islogical(x));
addParameter(p,'ReferenceFrequenciesHz',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'ReferenceReceiverHz',[],@(x) isempty(x)||(isnumeric(x)&&isscalar(x)));
addParameter(p,'ImagingReceiverHz',[],@(x) isempty(x)||(isnumeric(x)&&isscalar(x)));
addParameter(p,'DwellSeconds',[],@(x) isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'FrequencySign',1,@(x) isscalar(x)&&any(x==[1,-1]));
addParameter(p,'Verbose',false,@(x) islogical(x)||isscalar(x));
addParameter(p,'Plot',false,@(x) islogical(x)||isscalar(x));
addParameter(p,'ErrorOnSmallInput',false,@(x) islogical(x)||isscalar(x));
parse(p,varargin{:}); opt=p.Results;

D=opt.SpatialDimensions; nd=ndims(kspace); originalShape=size(kspace);
if isempty(opt.CoilDimension)
    if nd>D, cdim=nd; else, cdim=0; end
else
    cdim=opt.CoilDimension;
end
if cdim>nd || (cdim~=0 && cdim<=D)
    error('rs2d_remove_coherent_tones:CoilDimension','Invalid coil dimension.');
end
nx=size(kspace,1); ny=size(kspace,2);
if D==3, nz=size(kspace,3); else, nz=1; end
kspaceCorrected=kspace;
info=struct('version','5.0','status','not_processed','skipReason','', ...
    'spatialDimensions',D,'coilDimension',cdim,'frames',0, ...
    'candidateReadoutPositions',[],'refinedReadoutPositions',[], ...
    'candidateCoherences',[],'candidateThetaSpreads',[], ...
    'numAccepted',0,'frameSummary',[],'toneTests',[], ...
    'acquisitionIndex',[],'envelopeKnotTimes',[], ...
    'sliceBackgroundReduction',[],'sliceAnatomyChange',[], ...
    'ignoredOptions',p.Unmatched,'settings',opt);
if nd>=5 && size(kspace,3)>1 && ismember('SpatialDimensions',p.UsingDefaults)
    info.status='skipped_ambiguous_dimensions';
    info.skipReason='Specify SpatialDimensions for ambiguous 5D data.';return;
end
if D==3 && nd>=5 && ismember('CoilDimension',p.UsingDefaults)
    info.status='skipped_ambiguous_dimensions';
    info.skipReason='Specify CoilDimension (0 for no receiver dimension).';return;
end
if nx<16||ny<16||(D==3&&nz<2)
    if opt.ErrorOnSmallInput
        error('rs2d_remove_coherent_tones:SmallInput','Input is too small.');
    end
    info.status='skipped_small_input'; info.skipReason='Small/nonimaging array unchanged.';return;
end
if opt.MaxTones==0
    info.status='no_correction';info.skipReason='MaxTones=0.';return;
end
if any(~isfinite(real(kspace(:))))||any(~isfinite(imag(kspace(:))))
    error('rs2d_remove_coherent_tones:Nonfinite','Input contains nonfinite values.');
end
if cdim==0
    extra=(D+1):nd; perm=1:nd; nc=1;
else
    extra=setdiff((D+1):nd,cdim,'stable');
    perm=[1:D extra cdim];nc=size(kspace,cdim);
end
if isempty(extra), nf=1; else, nf=prod(originalShape(extra)); end
info.frames=nf;
K=reshape(double(permute(kspace,perm)),[nx ny nz nf nc]);
P=ny*nz;
if ~isempty(opt.AcquisitionIndex)
    A=reshape(double(opt.AcquisitionIndex),[],1);
    if numel(A)~=P||any(~isfinite(A))||any(A~=fix(A))||~isequal(sort(A),(0:P-1)')
        error('rs2d_remove_coherent_tones:AcquisitionIndex', ...
            'Expected a permutation of the integers 0..ny*nz-1.');
    end
    A=reshape(A,ny,nz);
else
    phaseOrder=char(opt.PhaseOrder);
    if strcmpi(phaseOrder,'auto')
        if D==3,phaseOrder='centric-alternating';else,phaseOrder='linear';end
    end
    A=makeAcquisitionIndex(ny,nz,phaseOrder,logical(opt.KZFast));
end
info.acquisitionIndex=A;
pidx=double(A(:));
ky=(0:ny-1)'; tail=max(4,floor(0.30*ny));
outer=(ky<tail)|(ky>=ny-tail);
out=repmat(outer,1,nz);out=out(:);
train=out & (mod(pidx,4)~=opt.ValidationFold);
val=out & (mod(pidx,4)==opt.ValidationFold);
if nnz(train)<16||nnz(val)<16
    info.status='skipped_small_input';info.skipReason='Insufficient peripheral readouts.';return;
end
sf=unique(max(1,min(nf,round(linspace(1,nf,min(nf,opt.CandidateFrames))))));
if isempty(opt.SupportMask)
    support=makeSupport(K,D,sf);
else
    support=logical(opt.SupportMask);
    if numel(support)~=nx*ny*nz
        error('rs2d_remove_coherent_tones:SupportMask','Image mask dimensions do not match spatial dimensions.');
    end
    support=reshape(support,[nx ny nz]);
end
if ~any(support(:)) || mean(support(:))>.78
    info.status='skipped_unreliable_mask';
    info.skipReason='Support empty or covers >78% of FOV.';return;
end
% Shared readout-frequency dictionary (candidate selection across frames).
if isempty(opt.ToneReadoutPositions)
    [candidates,coh,spread]=findCandidates(K,A,sf,train,opt);
else
    candidates=double(opt.ToneReadoutPositions(:).');
    if any(~isfinite(candidates))||any(candidates<1)||any(candidates>nx)
        error('rs2d_remove_coherent_tones:ToneReadoutPositions', ...
            'Positions must be 1-based and in [1,nx].');
    end
    coh=nan(size(candidates));spread=nan(size(candidates));
end
if ~isempty(opt.ReferenceFrequenciesHz)
    if isempty(opt.DwellSeconds)||isempty(opt.ReferenceReceiverHz)||isempty(opt.ImagingReceiverHz)
        error('rs2d_remove_coherent_tones:ReferenceMetadata', ...
            'ReferenceFrequenciesHz requires both receiver frequencies and DwellSeconds.');
    end
    freqs=double(opt.ReferenceFrequenciesHz(:).');
    df=freqs+opt.ReferenceReceiverHz-opt.ImagingReceiverHz;
    fs=1/opt.DwellSeconds;
    rpos=mod(floor(nx/2)+opt.FrequencySign*df*nx/fs,nx)+1;
    candidates=[candidates,rpos];
    coh=[coh,nan(size(rpos))];spread=[spread,nan(size(rpos))];
end
keep=false(size(candidates));saved=[];
for j=1:numel(candidates)
    if isempty(saved)||all(abs(saved-candidates(j))>1.8)
        keep(j)=true;saved(end+1)=candidates(j); %#ok<AGROW>
    end
end
candidates=candidates(keep);coh=coh(keep);spread=spread(keep);
candidates=candidates(1:min(numel(candidates),opt.MaxCandidates));
coh=coh(1:numel(candidates));spread=spread(1:numel(candidates));
info.candidateReadoutPositions=candidates;
info.candidateCoherences=coh;info.candidateThetaSpreads=spread;
if isempty(candidates)
    info.status='no_candidates';info.skipReason='No sufficiently coherent readout frequencies.';return;
end
refined=candidates;
if opt.RefineReadout
    for j=1:numel(candidates)
        refined(j)=refineReadout(K,candidates(j),pidx,train,sf);
    end
end
info.refinedReadoutPositions=refined;
% Smooth envelope basis in actual acquisition time.
nb=min(opt.EnvelopeKnots,P);
if nb==1
    B=ones(P,1);
else
    knotPositions=linspace(0,P-1,nb);
    B=interp1(knotPositions,eye(nb),pidx,'linear');
end
info.envelopeKnotTimes=linspace(0,P-1,nb);
Bt=B(train,:);
L=Bt'*Bt;s=max(trace(L)/nb,eps);
if nb>=3
    D2=diff(eye(nb),2,1);
else
    D2=zeros(0,nb);
end
ridgeMatrix=L + opt.EnvelopeRidge*s*(D2'*D2)+1e-5*s*eye(nb);
work=K;nAccepted=0;
empty=struct('frame',0,'tone',0,'readoutPosition',0, ...
    'phaseCyclesPerReadout',0,'trainingCoherence',0, ...
    'heldOutImprovement',0,'envelopeCV',0,'strength',0, ...
    'accepted',false,'rejectionReason','', ...
    'envelopeCoefficients',[]);
nt=numel(refined);tests=repmat(empty,nf*nt,1);r=0;
frameSummary=zeros(nf,4);
sliceBackground=zeros(nf,nz);
sliceAnatomy=zeros(nf,nz);
for f=1:nf
    k0=reshape(K(:,:,:,f,:),[nx ny nz nc]);
    k=reshape(work(:,:,:,f,:),[nx ny nz nc]);
    mag0=imageMagnitude(k0,D);
    den=max(mean(mag0(support)),eps);
    bg0=sum(abs(mag0(~support)).^2);
    sliceDen=ones(1,nz);
    for zidx=1:nz
        zMask=support(:,:,zidx); zMag=mag0(:,:,zidx);
        if any(zMask(:)), sliceDen(zidx)=max(mean(zMag(zMask)),eps);end
    end
    na=0;
    for j=1:nt
        x=refined(j);accepted=false;strength=0;
        theta=NaN;coher=NaN;improvement=NaN;envelopeCV=NaN;
        reason='';coefficients=[];
        if na>=opt.MaxTones
            reason='max_accepted_tones';
        else
            y=projectReadout(k,x); % P x nc
            theta=fitThetaTraining(y,pidx,train);
            ph=exp(2i*pi*theta*pidx);
            z=bsxfun(@times,y,conj(ph));
            coef=ridgeMatrix\(Bt'*z(train,:));
            coefficients=coef;
            envelope=B*coef;
            pred=bsxfun(@times,ph,envelope);
            before=sum(abs(y(val,:)).^2,'all');
            after=sum(abs(y(val,:)-pred(val,:)).^2,'all');
            improvement=1-after/max(before,eps);
            meanz=mean(z(train,:),1);
            coher=sum(abs(meanz).^2)/max(mean(sum(abs(y(train,:)).^2,2)),eps);
            envelopeCV=std(sqrt(sum(abs(envelope).^2,2))) / ...
                max(mean(sqrt(sum(abs(envelope).^2,2))),eps);
            if improvement<opt.MinValidationReduction
                reason='failed_holdout';
            elseif coher<opt.MinTemporalCoherence
                reason='insufficient_temporal_coherence';
            else
                bn=(0:nx-1)'-floor(nx/2);
                bx=exp(-2i*pi*(x-(floor(nx/2)+1))*bn/nx);
                model=bsxfun(@times,reshape(bx,[nx 1 1 1]), ...
                    reshape(pred,[1 ny nz nc]));
                mcand=imageMagnitude(k-model,D);
                [err,sliceErr]=anatomyErrors(mcand,mag0,support,den,sliceDen);
                if err<=opt.MaxTotalAnatomyChange && ...
                        max(sliceErr)<=opt.MaxPerSliceAnatomyChange
                    strength=1;
                else
                    lo=0;hi=1;
                    for it=1:8
                        mid=(lo+hi)/2;
                        mt=imageMagnitude(k-mid*model,D);
                        [change,localChanges]=anatomyErrors(mt,mag0,support,den,sliceDen);
                        if change>opt.MaxTotalAnatomyChange || ...
                                max(localChanges)>opt.MaxPerSliceAnatomyChange
                            hi=mid;
                        else
                            lo=mid;
                        end
                    end
                    strength=lo;
                end
                if strength>=opt.MinCorrectionStrength
                    k=k-strength*model;na=na+1;nAccepted=nAccepted+1;
                    accepted=true;reason='accepted';
                else
                    strength=0;reason='anatomy_change_limit';
                end
            end
        end
        r=r+1;
        tests(r)=struct('frame',f,'tone',j,'readoutPosition',x, ...
            'phaseCyclesPerReadout',theta,'trainingCoherence',coher, ...
            'heldOutImprovement',improvement,'envelopeCV',envelopeCV, ...
            'strength',strength,'accepted',accepted,'rejectionReason',reason, ...
            'envelopeCoefficients',coefficients);
    end
    mag1=imageMagnitude(k,D);
    bg=1-sum(abs(mag1(~support)).^2)/max(bg0,eps);
    anatomy=mean(abs(mag1(support)-mag0(support)))/den;
    frameSummary(f,:)=[f na bg anatomy];
    [~,sliceChange]=anatomyErrors(mag1,mag0,support,den,sliceDen);
    sliceAnatomy(f,:)=sliceChange;
    for zidx=1:nz
        zMask=support(:,:,zidx);
        za=mag0(:,:,zidx);zb=mag1(:,:,zidx);
        out0=sum(abs(za(~zMask)).^2);
        out1=sum(abs(zb(~zMask)).^2);
        sliceBackground(f,zidx)=1-out1/max(out0,eps);
    end
    work(:,:,:,f,:)=reshape(k,[nx ny nz 1 nc]);
    if opt.Verbose
        fprintf('v5 frame %d/%d: %d tones, bg %.2f%%, anatomy %.2f%%, peak slice %.2f%%\n', ...
            f,nf,na,100*bg,100*anatomy,100*max(sliceChange));
    end
end
info.numAccepted=nAccepted;info.toneTests=tests(1:r);
info.frameSummary=frameSummary;
info.sliceBackgroundReduction=sliceBackground;
info.sliceAnatomyChange=sliceAnatomy;
if nAccepted>0,info.status='corrected';else,info.status='no_correction';end
if cdim==0
    restored=reshape(work,originalShape);
else
    restored=ipermute(reshape(work,originalShape(perm)),perm);
end
kspaceCorrected=cast(restored,'like',kspace);
if opt.Plot
    z=ceil(nz/2);
    orig=imageMagnitude(reshape(K(:,:,:,1,:),[nx ny nz nc]),D);
    corr=imageMagnitude(reshape(work(:,:,:,1,:),[nx ny nz nc]),D);
    lim=localQuantile(orig(:),.995);
    figure('Name','v5 local coherent interference');
    subplot(1,3,1);imagesc(orig(:,:,z).');axis image off;caxis([0 lim]);colormap gray;title('Original');
    subplot(1,3,2);imagesc(corr(:,:,z).');axis image off;caxis([0 lim]);title('Corrected');
    subplot(1,3,3);imagesc((corr(:,:,z)-orig(:,:,z)).');axis image off;colorbar;title('Difference');
end
end

function [globalChange,sliceChange]=anatomyErrors(candidate,original,mask,den,sliceDen)
% Mean relative magnitude changes (global and separately by kz slice).
% This is a safeguard, not a proof that T1/R1 signal has been preserved.
[nx,ny,nz]=size(original); %#ok<ASGLU>
if any(mask(:))
    globalChange=mean(abs(candidate(mask)-original(mask)))/max(den,eps);
else
    globalChange=0;
end
sliceChange=zeros(1,nz);
for z=1:nz
    zz=mask(:,:,z);
    if any(zz(:))
        aa=candidate(:,:,z);bb=original(:,:,z);
        sliceChange(z)=mean(abs(aa(zz)-bb(zz)))/max(sliceDen(z),eps);
    end
end
end

function A=makeAcquisitionIndex(ny,nz,phaseOrder,kzFast)
switch lower(strtrim(phaseOrder))
    case 'linear',ord=1:ny;
    case {'centric','centric-alternating','zigzag','centre-out'}
        c=floor(ny/2)+1;ord=c;
        for d=1:ny
            if c+d<=ny,ord(end+1)=c+d;end %#ok<AGROW>
            if c-d>=1,ord(end+1)=c-d;end %#ok<AGROW>
        end
    otherwise, error('rs2d_remove_coherent_tones:PhaseOrder','Unknown phase order.');
end
rank=zeros(ny,1);rank(ord)=0:ny-1;
if kzFast
    A=repmat(rank,1,nz)*nz+repmat(0:nz-1,ny,1);
else
    A=repmat(rank,1,nz)+repmat((0:nz-1)*ny,ny,1);
end
end

function [cand,score,spread]=findCandidates(K,A,sf,train,opt)
[nx,ny,nz,~,nc]=size(K);P=ny*nz;
M=numel(sf);coh=zeros(nx,M);freq=zeros(nx,M);
ind=double(A(:))+1;
for jj=1:M
    k=reshape(K(:,:,:,sf(jj),:),[nx P nc]);
    xr=fftshift(ifft(ifftshift(k,1),[],1),1);
    temp=zeros(nx,P,nc);
    % Discover frequencies using TRAINING lines only (not held-out lines).
    kselect=find(train);
    temp(:,ind(kselect),:)=xr(:,kselect,:);
    spec=fft(temp,[],2);
    pow=sum(abs(spec).^2,3);
    [mx,where]=max(pow,[],2);
    den=P*sum(sum(abs(temp).^2,3),2);
    coh(:,jj)=mx./max(den,eps);
    freq(:,jj)=(where-1)/P;
end
med=median(coh,2);
fref=angle(mean(exp(2i*pi*freq),2))/(2*pi);
delta=mod(freq-fref+0.5,1)-0.5;
spreadAll=max(abs(delta),[],2);
possible=find(med>=opt.MinCandidateCoherence & spreadAll<=opt.MaxThetaSpread);
[~,ix]=sort(med(possible),'descend');possible=possible(ix);
cand=[];score=[];spread=[];
for jj=1:numel(possible)
    x=possible(jj);
    if ~isempty(cand)&&min(abs(cand-x))<2,continue;end
    cand(end+1)=x;score(end+1)=med(x);spread(end+1)=spreadAll(x); %#ok<AGROW>
    if numel(cand)>=opt.MaxCandidates,break;end
end
end

function best=refineReadout(K,x0,pidx,train,sf)
[nx,ny,nz,~,nc]=size(K);
xx=x0+(-.70:.14:.70);scores=zeros(size(xx));
sf=sf(1:min(3,numel(sf)));
for j=1:numel(xx)
    for f=sf
        k=reshape(K(:,:,:,f,:),[nx ny nz nc]);
        y=projectReadout(k,xx(j));
        theta=fitThetaTraining(y,pidx,train);
        ph=exp(-2i*pi*theta*pidx(train));
        b=mean(bsxfun(@times,y(train,:),ph),1);
        scores(j)=scores(j)+sum(abs(b).^2);
    end
end
[~,ii]=max(scores);best=xx(ii);
if ii>1&&ii<numel(xx)
    a=scores(ii-1);b=scores(ii);c=scores(ii+1);
    den=a-2*b+c;
    if abs(den)>eps&&den<0
        best=best+0.14*max(-1,min(1,0.5*(a-c)/den));
    end
end
best=max(1,min(nx,best));
end

function y=projectReadout(k,x)
[nx,ny,nz,nc]=size(k);P=ny*nz;
n=(0:nx-1)'-floor(nx/2);
b=exp(-2i*pi*(x-(floor(nx/2)+1))*n/nx);
y=(conj(b).'*reshape(k,nx,P*nc))/nx;
y=reshape(y,[P nc]);
end

function theta=fitThetaTraining(y,pidx,train)
[P,nc]=size(y);nfft=16384;
q=pidx(train)+1;
temp=zeros(P,nc);
temp(q,:)=y(train,:);
f=fft(temp,nfft,1);
spec=sum(abs(f).^2,2);
[~,where]=max(spec);
theta0=(where-1)/nfft;
times=pidx(train);
subset=y(train,:);
obj=@(theta) -sum(abs(sum(bsxfun(@times,subset,exp(-2i*pi*theta*times)),1)).^2);
theta=fminbnd(obj,theta0-1.0/nfft,theta0+1.0/nfft, ...
    optimset('TolX',5e-9,'Display','off'));
theta=mod(theta,1);
end

function m=imageMagnitude(K,D)
[nx,ny,nz,nc]=size(K);I=K;
for d=1:D,I=ifftshift(I,d);end
for d=1:D,I=ifft(I,[],d);end
for d=1:D,I=fftshift(I,d);end
m=sqrt(sum(abs(I).^2,4));m=reshape(m,[nx ny nz]);
end

function support=makeSupport(K,D,sf)
[nx,ny,nz,~,nc]=size(K);
a=zeros(nx,ny,nz);
for f=sf
    k=reshape(K(:,:,:,f,:),[nx ny nz nc]);
    m=imageMagnitude(k,D);
    a=a+m/max(localQuantile(m(:),.995),eps)/numel(sf);
end
a=movmedian(a,11,1);a=movmedian(a,5,2);
a=gaussianSmooth(a,1,2.5);a=gaussianSmooth(a,2,2);
if D==3,a=gaussianSmooth(a,3,.6);end
support=a>.18*localQuantile(a(:),.995);
if D==3
    support=convn(double(support),ones(3,3,3),'same')>0;
else
    support=convn(double(support),ones(3,3),'same')>0;
end
end

function a=gaussianSmooth(a,dim,sigma)
r=ceil(3*sigma);u=-r:r;g=exp(-.5*(u/sigma).^2);g=g/sum(g);
sz=ones(1,max(3,ndims(a)));sz(dim)=numel(g);
a=convn(a,reshape(g,sz),'same');
end

function q=localQuantile(v,p)
v=sort(double(v(:)));j=1+(numel(v)-1)*p;
a=floor(j);b=ceil(j);
q=v(a)+(j-a)*(v(b)-v(a));
end
