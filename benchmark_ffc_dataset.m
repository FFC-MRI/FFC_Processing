function report = benchmark_ffc_dataset(scans)
%BENCHMARK_FFC_DATASET Measure preload and repeated scan selection numerics.
% scans = open_mrifile(); report = benchmark_ffc_dataset(scans);
% Measures numerical reuse, not UI event latency. Does not mutate input scans.
    ffc_setup();
    options = struct('reconstruct',true,'reprocess',false);
    started = tic;
    [ready,~,failures] = ffc_preprocess_scans(scans,options);
    report = struct('loadSeconds',toc(started),'sequenceCount',numel(ready), ...
        'failedSequenceCount',numel(failures),'matlabVersion',version);
    if isempty(ready), error('benchmark_ffc_dataset:NoScans','No scans reconstructed.'); end
    generations = cellfun(@(s) s.reconstruction_generation,ready);
    times = zeros(1,3*numel(ready));
    order = repmat(1:numel(ready),1,3);
    for k = 1:numel(order)
        index = order(k); started = tic;
        ready{index} = ffc_normalize_viewer_scan(ready{index});
        ready{index} = buildimages(ready{index},false);
        times(k) = toc(started);
    end
    assert(isequal(generations,cellfun(@(s) s.reconstruction_generation,ready)), ...
        'Changing selection unexpectedly rebuilt an image.');
    report.medianSelectionSeconds = median(times);
    report.maxSelectionSeconds = max(times);
    report.unexpectedRebuilds = 0;
    report.note = 'Numerical selection/cache checks only; loading includes any initial reconstruction/cache writes.';
    disp(report);
end
