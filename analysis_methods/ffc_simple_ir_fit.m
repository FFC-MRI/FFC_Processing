function [result, R1, rsquare] = ffc_simple_ir_fit(timesMs, signedSignal)
%FFC_SIMPLE_IR_FIT Ordinary unweighted offset-plus-exponential IR fit.
%
% The fitted model is
%       y(t) = A*exp(-t/T1) + C
% which is algebraically identical to
%       y(t) = (1-2*m0)*exp(-t/t1) + c.
%
% For a trial T1, A and C are obtained by ordinary linear least squares.
% Only the positive scalar T1 is searched numerically. This gives the same
% least-squares objective as a three-parameter nonlinear fit without adding
% phase, noise-floor, robust-weighting or joint-field parameters.

    t = double(timesMs(:));
    y = double(signedSignal(:));
    valid = isfinite(t) & isfinite(y);
    t = t(valid);
    y = y(valid);
    if numel(t) < 3
        error('ffc_simple_ir_fit:TooFewPoints', ...
            'At least three finite time points are required.');
    end

    [t, order] = sort(t);
    y = y(order);
    if numel(unique(t)) < 2
        error('ffc_simple_ir_fit:NoTimeVariation', ...
            'At least two distinct evolution times are required.');
    end

    timeScale = max(abs(t));
    if ~isfinite(timeScale) || timeScale <= 0
        timeScale = 1;
    end
    lowerT1 = max(1e-6, eps(timeScale));
    upperT1 = max(1e6, 1e4 .* timeScale);
    persistent options
    if isempty(options), options = optimset('Display', 'off', 'TolX', 1e-9); end
    objective = @(logT1) linearFitAtT1(exp(logT1), t, y);
    [bestLogT1, bestSse] = fminbnd(objective, ...
        log(lowerT1), log(upperT1), options);
    T1ms = exp(bestLogT1);
    [~, coefficients, prediction] = linearFitAtT1(T1ms, t, y);

    amplitude = coefficients(1);
    offset = coefficients(2);
    residual = y - prediction;
    sst = sum((y - mean(y)).^2);
    if sst > 0
        rsquare = 1 - bestSse ./ sst;
    else
        rsquare = NaN;
    end
    dfe = max(numel(y) - 3, 0);

    R1 = 1000 ./ T1ms;
    if nargout > 1
        result = T1ms;
        return;
    end
    result = struct();
    result.modelEquation = 'A*exp(-t/T1)+C';
    result.timesMs = t;
    result.signedSignal = y;
    result.prediction = prediction;
    result.residual = residual;
    result.A = amplitude;
    result.C = offset;
    result.m0 = (1 - amplitude) ./ 2;
    result.T1ms = T1ms;
    result.R1 = 1000 ./ T1ms;
    result.sse = bestSse;
    result.rsquare = rsquare;
    result.rmse = sqrt(bestSse ./ numel(y));
    result.dfe = dfe;
end


function [sse, coefficients, prediction] = linearFitAtT1(T1ms, t, y)
    if ~isfinite(T1ms) || T1ms <= 0
        sse = realmax('double');
        coefficients = [NaN; NaN];
        prediction = nan(size(y));
        return;
    end
    design = [exp(-t ./ T1ms), ones(size(t))];
    coefficients = design \ y;
    prediction = design * coefficients;
    residual = y - prediction;
    sse = sum(residual.^2);
    if ~isfinite(sse)
        sse = realmax('double');
    end
end
