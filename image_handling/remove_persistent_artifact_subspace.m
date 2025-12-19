function [Kcorr, info] = remove_persistent_artifact_subspace(K, opts)
% Remove narrowband periodic RF pickup along readout (dim 1),
% but only apply subtraction strongly in outer kx with smooth taper to zero at centre.
%
% Assumes: dim1=readout, dim2=phase, last dim=coils. Other dims allowed.

arguments
    K
    opts.guardBins (1,1) double = 12          % exclude near-DC bins for frequency estimate
    opts.learnOuterFrac (1,1) double = 0.30   % samples at each edge used to estimate tone
    opts.applyOuterFrac (1,1) double = 0.45   % region where subtraction is applied (tapered)
    opts.maxLinesForFreq (1,1) double = 4096
    opts.forceBin (1,1) double = NaN
end

K = complex(K);
sz = size(K);
Nro = sz(1);
Nc  = sz(end);
other = prod(sz(2:end-1));      % collapses ky and any other dims except coils

Krs = reshape(K, Nro, other, Nc);

% --- masks along readout ---
nLearn = max(8, round(opts.learnOuterFrac * Nro));
learnMask = false(Nro,1);
learnMask(1:nLearn) = true;
learnMask(end-nLearn+1:end) = true;

% Tapered application window: 1 at edges -> 0 at centre
nApply = max(8, round(opts.applyOuterFrac * Nro));
w = zeros(Nro,1);
w(1:nApply) = 1;
w(end-nApply+1:end) = 1;
% raised-cosine taper from edge region into centre
% make w smooth so we don't introduce ringing/black bars
taperLen = max(8, round(0.15 * Nro));
mid1 = nApply - taperLen + 1;
mid2 = nApply;
if mid1 > 1
    t = (0:taperLen-1).'/(taperLen-1);
    win = 0.5*(1+cos(pi*t));     % 1 -> 0
    w(mid1:mid2) = win;
    w(end-mid2+1:end-mid1+1) = flipud(win);
end

% --- estimate tone bin along readout FFT using only outer samples (learnMask) ---
if isnan(opts.forceBin)
    nUse = min(other, opts.maxLinesForFreq);
    idx = round(linspace(1, other, nUse));

    Pacc = zeros(Nro,1);
    for c = 1:Nc
        X = Krs(:, idx, c);
        X(~learnMask,:) = 0;                 % suppress anatomy-heavy centre
        Sf = fft(X, [], 1);
        Pacc = Pacc + mean(abs(Sf), 2);
    end

    g = min(opts.guardBins, floor(Nro/4));
    Pacc(1:g) = 0; Pacc(end-g+1:end) = 0;
    [~, fbin] = max(Pacc);
else
    fbin = round(opts.forceBin);
end

% normalized frequency cycles/sample in [-0.5,0.5)
k = fbin - 1;
f = k / Nro;
if f >= 0.5, f = f - 1; end

% basis
n = (0:Nro-1).';
c1 = cos(2*pi*f*n);
s1 = sin(2*pi*f*n);

A = [c1(learnMask), s1(learnMask)];
AtA_inv = inv(A' * A);          % 2x2

Kcorr_rs = Krs;

for c = 1:Nc
    X = Krs(:,:,c);             % [Nro x other]
    Xm = X(learnMask,:);        % learn region only

    % fit real & imag separately
    ar = AtA_inv * (A' * real(Xm));
    ai = AtA_inv * (A' * imag(Xm));

    tone = complex(c1*ar(1,:) + s1*ar(2,:), c1*ai(1,:) + s1*ai(2,:) );

    % APPLY SUBTRACTION WITH TAPER (critical!)
    Kcorr_rs(:,:,c) = X - (w .* tone);
end

Kcorr = reshape(Kcorr_rs, sz);

info = struct();
info.f_bin = fbin;
info.f_cycles_per_sample = f;
info.learnOuterFrac = opts.learnOuterFrac;
info.applyOuterFrac = opts.applyOuterFrac;
info.note = "Tapered tone subtraction along readout to avoid creating y-localised black bars.";
end
