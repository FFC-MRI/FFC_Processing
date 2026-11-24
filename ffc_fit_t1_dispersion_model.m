function fitResult = ffc_fit_t1_dispersion_model(fields_mT, R1, modelName, quadrupoleMode, options)
%FFC_FIT_T1_DISPERSION_MODEL Fit an R1 NMRD/FFC dispersion model.
%
%   fitResult = ffc_fit_t1_dispersion_model(fields_mT,R1,modelName,...
%       quadrupoleMode,options)
%
% The fit is always performed in R1 because independent relaxation
% contributions are additive in relaxation rate. The returned structure also
% contains T1 predictions for display.
%
% Supported background models:
%   'Power law (biological)'          Rinf + A*f^(-beta)
%   'BPP single tau (standard)'       Rinf + DeltaR*g(f,tau)
%   'BPP two tau (biological)'        Rinf + D1*g(f,tau1)+D2*g(f,tau2)
%   'Power law + BPP (biological)'    Rinf + A*f^(-beta)+D*g(f,tau)
%
% where f is the proton Larmor frequency in MHz and
% g = [J(w)+4J(2w)]/(5*tau), J(w)=tau/(1+(w*tau)^2), so g(0)=1.
%
% Optional empirical 14N quadrupolar-enhanced-relaxation (QRE) terms are
% added as positive Lorentzians in R1. Their resonance centres are fixed to
% literature values by default because typical FCI imaging dispersion curves
% are too sparsely sampled to robustly determine centre frequency as well as
% amplitude. The shared HWHM can optionally be fitted.
%
% options fields:
%   qHWHMMHz       fixed/start shared Lorentzian HWHM (default 0.20 MHz)
%   fitQWidth      fit a shared QRE HWHM (default false)
%   nStarts        deterministic multistart count (default 24)
%
% This routine uses LSQCURVEFIT when available and otherwise falls back to a
% bounded parameter transform with FMINSEARCH.

    if nargin < 5 || isempty(options), options = struct(); end
    if nargin < 4 || isempty(quadrupoleMode), quadrupoleMode = 'None'; end
    if nargin < 3 || isempty(modelName), modelName = 'Power law (biological)'; end

    fields_mT = double(fields_mT(:));
    R1 = double(R1(:));
    valid = isfinite(fields_mT) & fields_mT > 0 & isfinite(R1) & R1 > 0;
    fields_mT = fields_mT(valid);
    R1 = R1(valid);
    if numel(R1) < 4
        error('ffc_fit_t1_dispersion_model:TooFewPoints', ...
            'At least four valid positive R1 values are required.');
    end

    gammaMHzPerT = 42.57747892;
    frequencyMHz = gammaMHzPerT .* fields_mT ./ 1000;

    qHWHM = optionNumeric(options,'qHWHMMHz',0.20);
    qHWHM = min(max(qHWHM,0.03),1.0);
    fitQWidth = optionLogical(options,'fitQWidth',false);
    nStarts = max(1,round(optionNumeric(options,'nStarts',24)));

    [base, p0, lb, ub, paramNames, paramUnits] = ...
        baseModelDefinition(modelName,frequencyMHz,R1);
    [qCenters, qLabel] = quadrupoleDefinition(quadrupoleMode);
    qCount = numel(qCenters);

    yScale = max(max(R1),1);
    if qCount > 0
        qAmp0 = repmat(max(rangeFinite(R1)/max(qCount,1),0.05*yScale),1,qCount);
        qAmpLb = zeros(1,qCount);
        qAmpUb = repmat(max(5*yScale,5),1,qCount);
        p0 = [p0 qAmp0]; %#ok<AGROW>
        lb = [lb qAmpLb]; %#ok<AGROW>
        ub = [ub qAmpUb]; %#ok<AGROW>
        for k=1:qCount
            paramNames{end+1} = sprintf('Q%d_amplitude',k); %#ok<AGROW>
            paramUnits{end+1} = 's^-1'; %#ok<AGROW>
        end
        if fitQWidth
            p0(end+1) = qHWHM;
            lb(end+1) = 0.03;
            ub(end+1) = 0.80;
            paramNames{end+1} = 'Q_shared_HWHM';
            paramUnits{end+1} = 'MHz';
        end
    end

    modelFun = @(p,f) evaluateFullModel(p,f,base,qCenters,qHWHM,fitQWidth);
    nParams = numel(p0);
    nPoints = numel(R1);
    if nPoints <= nParams
        error('ffc_fit_t1_dispersion_model:Overparameterised', ...
            ['This model has %d fitted parameters but only %d valid field ', ...
             'points. Select a simpler background model or fewer quadrupole terms.'], ...
            nParams,nPoints);
    end

    [bestP, bestSSE, exitflag, solverName] = boundedMultistart( ...
        modelFun,frequencyMHz,R1,p0,lb,ub,nStarts);

    prediction = modelFun(bestP,frequencyMHz);
    residual = R1-prediction;
    sse = sum(residual.^2);
    rmse = sqrt(mean(residual.^2));
    sst = sum((R1-mean(R1)).^2);
    if sst > 0, rsquare = 1-sse/sst; else, rsquare = NaN; end
    if nPoints > nParams+1
        adjustedR2 = 1-(1-rsquare)*(nPoints-1)/(nPoints-nParams-1);
    else
        adjustedR2 = NaN;
    end
    mseForIC = max(sse/nPoints,realmin);
    aic = nPoints*log(mseForIC)+2*nParams;
    if nPoints > nParams+1
        aicc = aic + (2*nParams*(nParams+1))/(nPoints-nParams-1);
    else
        aicc = NaN;
    end
    bic = nPoints*log(mseForIC)+nParams*log(nPoints);

    denseFrequencyMHz = logspace(log10(min(frequencyMHz)*0.85), ...
        log10(max(frequencyMHz)*1.15),500).';
    denseR1 = modelFun(bestP,denseFrequencyMHz);
    denseT1ms = 1000./denseR1;

    qWidthUsed = qHWHM;
    if qCount > 0 && fitQWidth
        qWidthUsed = bestP(end);
    end

    qStart = base.nParameters + 1;
    if qCount > 0
        qAmplitudes = bestP(qStart:(qStart+qCount-1));
        qAreas = pi .* qAmplitudes .* qWidthUsed;
    else
        qAmplitudes = [];
        qAreas = [];
    end

    % Reorder the two-tau BPP model for consistent reporting.
    if strcmp(base.key,'bpp2')
        % [Rinf,D1,tau1,D2,tau2]
        if bestP(3) > bestP(5)
            bestP([2 4]) = bestP([4 2]);
            bestP([3 5]) = bestP([5 3]);
            prediction = modelFun(bestP,frequencyMHz);
            residual = R1-prediction;
            denseR1 = modelFun(bestP,denseFrequencyMHz);
            denseT1ms = 1000./denseR1;
        end
    end

    fitResult = struct();
    fitResult.modelName = char(modelName);
    fitResult.modelKey = base.key;
    fitResult.modelEquation = base.equation;
    fitResult.fitDomain = 'R1';
    fitResult.quadrupoleMode = char(quadrupoleMode);
    fitResult.quadrupoleLabel = qLabel;
    fitResult.quadrupoleCentersMHz = qCenters(:).';
    fitResult.quadrupoleHWHMMHz = qWidthUsed;
    fitResult.quadrupoleWidthFitted = logical(fitQWidth && qCount>0);
    fitResult.quadrupoleAmplitudes_s_1 = qAmplitudes(:).';
    fitResult.quadrupoleAreas_s_1_MHz = qAreas(:).';
    fitResult.totalQuadrupoleArea_s_1_MHz = sum(qAreas);
    fitResult.parameters = bestP(:).';
    fitResult.parameterNames = paramNames;
    fitResult.parameterUnits = paramUnits;
    fitResult.fields_mT = fields_mT(:).';
    fitResult.frequencyMHz = frequencyMHz(:).';
    fitResult.observedR1 = R1(:).';
    fitResult.observedT1ms = (1000./R1(:)).';
    fitResult.predictedR1 = prediction(:).';
    fitResult.predictedT1ms = (1000./prediction(:)).';
    fitResult.residualR1 = residual(:).';
    fitResult.denseFrequencyMHz = denseFrequencyMHz(:).';
    fitResult.denseField_mT = denseFrequencyMHz(:).' ./ gammaMHzPerT .* 1000;
    fitResult.denseR1 = denseR1(:).';
    fitResult.denseT1ms = denseT1ms(:).';
    fitResult.sse = sse;
    fitResult.rmse = rmse;
    fitResult.rsquare = rsquare;
    fitResult.adjustedRsquare = adjustedR2;
    fitResult.AIC = aic;
    fitResult.AICc = aicc;
    fitResult.BIC = bic;
    fitResult.nPoints = nPoints;
    fitResult.nParameters = nParams;
    fitResult.degreesOfFreedom = nPoints-nParams;
    fitResult.exitflag = exitflag;
    fitResult.solver = solverName;
    fitResult.gammaMHzPerT = gammaMHzPerT;
    fitResult.notes = ['Quadrupole terms are empirical positive Lorentzians ', ...
        'in R1; fixed centres do not constitute a full microscopic QRE model.'];
end

function [base,p0,lb,ub,names,units] = baseModelDefinition(modelName,f,R1)
    name = lower(strtrim(char(modelName)));
    yMin = min(R1);
    yMax = max(R1);
    yRange = max(yMax-yMin,0.1*yMax);
    yCap = max(10*yMax,10);
    fMed = median(f);
    tauGuessUs = 1/(2*pi*max(fMed,1e-6)); % because f is MHz -> tau in us
    tauGuessUs = min(max(tauGuessUs,0.005),500);

    if contains(name,'power law + bpp')
        base.key='powerbpp';
        base.nParameters=5;
        base.equation='R1=Rinf+A*(f/1MHz)^(-beta)+DeltaR*g_BPP(f,tau)';
        base.fun=@(p,x) p(1)+p(2).*max(x,1e-9).^(-p(3)) + ...
            p(4).*bppShape(x,p(5));
        p0=[max(0,0.5*yMin), max(yRange,0.1), 0.35, 0.5*yRange, tauGuessUs];
        lb=[0,0,0,0,1e-4];
        ub=[yCap,10*yCap,2.0,10*yCap,1e4];
        names={'Rinf','A_power','beta','DeltaR_BPP','tau_BPP'};
        units={'s^-1','s^-1 at 1 MHz','','s^-1','us'};
    elseif contains(name,'two tau') || contains(name,'two-tau')
        base.key='bpp2';
        base.nParameters=5;
        base.equation='R1=Rinf+D1*g_BPP(f,tau1)+D2*g_BPP(f,tau2)';
        base.fun=@(p,x) p(1)+p(2).*bppShape(x,p(3)) + ...
            p(4).*bppShape(x,p(5));
        p0=[max(0,0.7*yMin),0.55*yRange,max(tauGuessUs/4,1e-4), ...
            0.45*yRange,min(tauGuessUs*4,1e4)];
        lb=[0,0,1e-4,0,1e-4];
        ub=[yCap,10*yCap,1e4,10*yCap,1e4];
        names={'Rinf','DeltaR_fast','tau_fast','DeltaR_slow','tau_slow'};
        units={'s^-1','s^-1','us','s^-1','us'};
    elseif contains(name,'bpp') || contains(name,'single tau') || contains(name,'single-tau')
        base.key='bpp1';
        base.nParameters=3;
        base.equation='R1=Rinf+DeltaR*g_BPP(f,tau)';
        base.fun=@(p,x) p(1)+p(2).*bppShape(x,p(3));
        p0=[max(0,0.8*yMin), max(yRange,0.1), tauGuessUs];
        lb=[0,0,1e-4];
        ub=[yCap,10*yCap,1e4];
        names={'Rinf','DeltaR','tau'};
        units={'s^-1','s^-1','us'};
    else
        base.key='power';
        base.nParameters=3;
        base.equation='R1=Rinf+A*(f/1MHz)^(-beta)';
        base.fun=@(p,x) p(1)+p(2).*max(x,1e-9).^(-p(3));
        p0=[max(0,0.5*yMin), max(yRange,0.1), 0.35];
        lb=[0,0,0];
        ub=[yCap,10*yCap,2.0];
        names={'Rinf','A_power','beta'};
        units={'s^-1','s^-1 at 1 MHz',''};
    end
end

function [centres,label] = quadrupoleDefinition(mode)
    textMode = lower(strtrim(char(mode)));
    if contains(textMode,'three') || contains(textMode,'all')
        centres=[0.70 2.08 2.76];
        label='14N three-peak empirical QRE';
    elseif contains(textMode,'high') || contains(textMode,'pair') || contains(textMode,'2.1')
        centres=[2.08 2.76];
        label='14N high-field pair';
    elseif contains(textMode,'broad')
        centres=2.0;
        label='Broad approximately 2 MHz QRE';
    else
        centres=[];
        label='None';
    end
end

function y = evaluateFullModel(p,f,base,qCenters,qHWHM,fitQWidth)
    y = base.fun(p(1:base.nParameters),f);
    if isempty(qCenters), return; end
    qCount=numel(qCenters);
    amp = p((base.nParameters+1):(base.nParameters+qCount));
    width=qHWHM;
    if fitQWidth
        width=p(base.nParameters+qCount+1);
    end
    for k=1:qCount
        y = y + amp(k)./(1+((f-qCenters(k))./width).^2);
    end
end

function g = bppShape(fMHz,tauUs)
    tau = max(tauUs,realmin).*1e-6;
    w = 2*pi.*max(fMHz,0).*1e6;
    j1 = tau./(1+(w.*tau).^2);
    j2 = tau./(1+(2*w.*tau).^2);
    g = (j1+4*j2)./(5*tau);
end

function [bestP,bestSSE,exitflag,solverName] = boundedMultistart(fun,x,y,p0,lb,ub,nStarts)
    p0=min(max(p0,lb+1e-10.*max(1,abs(lb))),ub-1e-10.*max(1,abs(ub)));
    bestP=p0; bestSSE=Inf; exitflag=NaN; solverName='fminsearch bounded transform';

    starts=deterministicStarts(p0,lb,ub,nStarts);
    if exist('lsqcurvefit','file')==2
        solverName='lsqcurvefit multistart';
        try
            opts=optimoptions('lsqcurvefit','Display','off', ...
                'MaxIterations',2000,'MaxFunctionEvaluations',20000, ...
                'FunctionTolerance',1e-11,'StepTolerance',1e-11);
        catch
            opts=optimset('Display','off','MaxIter',2000,'MaxFunEvals',20000, ...
                'TolFun',1e-11,'TolX',1e-11);
        end
        for s=1:size(starts,1)
            try
                [p,resnorm,~,ef]=lsqcurvefit(fun,starts(s,:),x,y,lb,ub,opts);
                if isfinite(resnorm) && resnorm<bestSSE
                    bestSSE=resnorm; bestP=p; exitflag=ef;
                end
            catch
            end
        end
    end

    if ~isfinite(bestSSE)
        for s=1:size(starts,1)
            z0=toUnbounded(starts(s,:),lb,ub);
            objective=@(z) sum((y-fun(fromUnbounded(z,lb,ub),x)).^2);
            try
                opts=optimset('Display','off','MaxIter',5000,'MaxFunEvals',30000, ...
                    'TolFun',1e-11,'TolX',1e-9);
                [z,val,ef]=fminsearch(objective,z0,opts);
                p=fromUnbounded(z,lb,ub);
                if isfinite(val) && val<bestSSE
                    bestSSE=val; bestP=p; exitflag=ef;
                end
            catch
            end
        end
    end

    if ~isfinite(bestSSE)
        error('ffc_fit_t1_dispersion_model:FitFailed', ...
            'The dispersion model could not be fitted.');
    end
end

function starts=deterministicStarts(p0,lb,ub,nStarts)
    n=numel(p0);
    starts=zeros(nStarts,n);
    starts(1,:)=p0;
    if nStarts==1, return; end
    % Golden-ratio sequence gives repeatable coverage without modifying RNG.
    phi=(sqrt(5)-1)/2;
    for s=2:nStarts
        frac=mod((1:n).*phi*s,1);
        candidate=zeros(1,n);
        for j=1:n
            if lb(j)>0 && ub(j)/lb(j)>100
                % Correlation times span many decades: cover them
                % logarithmically rather than placing almost every start at
                % the slow end of the interval.
                candidate(j)=exp(log(lb(j))+frac(j)*(log(ub(j))-log(lb(j))));
            else
                candidate(j)=lb(j)+frac(j)*(ub(j)-lb(j));
            end
        end
        % Keep half the starts closer to the physically informed initial guess.
        if mod(s,2)==0
            candidate=0.65*p0+0.35*candidate;
        end
        starts(s,:)=min(max(candidate,lb),ub);
    end
end

function z=toUnbounded(p,lb,ub)
    q=(p-lb)./(ub-lb);
    q=min(max(q,1e-8),1-1e-8);
    z=log(q./(1-q));
end

function p=fromUnbounded(z,lb,ub)
    q=1./(1+exp(-max(min(z,40),-40)));
    p=lb+(ub-lb).*q;
end

function v=optionNumeric(options,name,defaultValue)
    v=defaultValue;
    if isstruct(options) && isfield(options,name) && ...
            isnumeric(options.(name)) && isscalar(options.(name)) && ...
            isfinite(options.(name))
        v=double(options.(name));
    end
end

function v=optionLogical(options,name,defaultValue)
    v=defaultValue;
    if isstruct(options) && isfield(options,name) && isscalar(options.(name))
        try, v=logical(options.(name)); catch, end
    end
end

function r=rangeFinite(x)
    x=x(isfinite(x));
    if isempty(x), r=0; else, r=max(x)-min(x); end
end
