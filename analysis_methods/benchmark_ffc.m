function report = benchmark_ffc(scan, outputFile)
%BENCHMARK_FFC Measure actual reconstruction and reuse on a representative scan.
%   report = benchmark_ffc(scan,'my_ffc_timings.mat')
% Use a preprocessed ImageReconCore, or supply one straight from open_mrifile.
% Does not change the input scan (ImageReconCore is a value class).
    if nargin < 2, outputFile = ''; end
    ffc_setup();
    if isempty(scan.originalcomplexkspace), scan = preprocessing(scan); end
    report = struct('matlabVersion',version,'computer',computer, ...
        'date',datestr(now,30),'sourceSize',size(scan.originalcomplexkspace), ...
        'fftSize',[scan.fft_size scan.fft_size_3d], ...
        'combination',scan.combination,'noiseWhitening',scan.noise_whitening, ...
        'selectedReceivers',scan.multichannel_recon);
    t = tic; ready = buildimages(scan); report.firstRebuildSeconds = toc(t);
    rebuild = zeros(1,3);
    for trial = 1:3
        t = tic; ready = buildimages(ready); rebuild(trial) = toc(t);
    end
    report.rebuildSeconds = rebuild;
    report.medianRebuildSeconds = median(rebuild);
    reuse = zeros(1,20); generation = ready.reconstruction_generation;
    for trial = 1:numel(reuse)
        t = tic; ready = buildimages(ready,false); reuse(trial) = toc(t);
    end
    assert(ready.reconstruction_generation == generation,'Cache unexpectedly rebuilt.');
    report.medianReuseSeconds = median(reuse);
    report.outputSize = size(ready.magimage);
    if isempty(ready.is_svs) || ~ready.is_svs
        planeTimes = zeros(1,30);
        volume = ready.magimage(:,:,:,1,1);
        for trial = 1:numel(planeTimes)
            t = tic;
            ffc_extract_anatomical_plane(volume,ready,'Axial',1,eye(3));
            planeTimes(trial) = toc(t);
        end
        report.medianPlaneSeconds = median(planeTimes);
    end
    report.note = ['Wall times on this MATLAB session; first rebuild includes JIT. ', ...
        'Three later rebuilds force all reconstruction work. Reuse is a distinct operation. ', ...
        'Display event latency and peak memory require separate measurement.'];
    disp(report);
    if ~isempty(outputFile), save(outputFile,'report'); end
end
