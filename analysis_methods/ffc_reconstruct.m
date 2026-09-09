function scan = ffc_reconstruct(scan, force)
%FFC_RECONSTRUCT Reconstruct one time/field volume at a time.
% 3D partitions stay together. Neither time, field nor coil is an FFT axis.
    if isempty(scan.originalcomplexkspace)
        error('ImageReconCore:NoPreprocessedKspace', ...
            'Run preprocessing successfully before buildimages.');
    end
    scan.complexkspace = scan.originalcomplexkspace;
    scan.fov_offsets_applied = false;
    scan = apply_fov_offset_phase_ramp(scan);
    source = scan.complexkspace;
    sz = size(source); sz(end+1:6) = 1;
    nc = sz(6); sel = scan.multichannel_recon;
    if isempty(sel)
        sel = 1:nc;
    elseif islogical(sel)
        if numel(sel) ~= nc
            error('buildimages:coilMaskWrongLength', 'Receiver mask must have %d elements.', nc);
        end
        sel = find(sel(:).');
    elseif isnumeric(sel)
        sel = unique(round(sel(:).'), 'stable');
        sel = sel(isfinite(sel));
    else
        error('buildimages:coilSelNotNumeric', 'Receiver selection must be numeric or logical.');
    end
    if isempty(sel), error('buildimages:coilSelEmptyAfterParse', 'Select at least one receiver.'); end
    if any(sel < 1 | sel > nc)
        error('buildimages:coilIndexOutOfRange', 'Receiver indices must lie in 1..%d.', nc);
    end
    scan.selected_coils_original = sel;
    nSel = numel(sel);
    whiten = ~isempty(scan.noise_whitening) && isscalar(scan.noise_whitening) && ...
        logical(scan.noise_whitening) && nc > 1;
    % Calibration depends on the corrected source, not FFT size or selection.
    calibrationKey = struct('source',scan.source_generation,'size',sz, ...
        'offsets',scan.fov_offset_info);
    if whiten
        if ~force && ~isempty(scan.whitening_cache) && ...
                isequaln(scan.whitening_cache_key, calibrationKey)
            info = scan.whitening_cache;
        else
            info = ffc_noise_calibration(source);
            scan.whitening_cache = info;
            scan.whitening_cache_key = calibrationKey;
        end
    else
        info = struct('method','disabled','nSamples',0);
        if nc == 1, info.method = 'singleAcquiredChannel'; end
    end
    info.nInputChannels = nc; info.nOutputChannels = nSel;
    info.selectedIndices = sel; info.selectedReceiverIds = sel - 1;
    scan.whitening_info = info;

    targetRead = round(double(scan.fft_size));
    if ~isscalar(targetRead) || ~isfinite(targetRead) || targetRead < 1
        error('buildimages:invalidFftSize', 'fft_size must be a positive finite scalar.');
    end
    target = [targetRead round(sz(2)*targetRead/sz(1)) sz(3)];
    if scan.TwoDimensional ~= 1 && ~isempty(scan.fft_size_3d)
        target(3) = round(double(scan.fft_size_3d));
    end
    if target(1) < sz(1), error('buildimages:fftSizeTooSmall', 'Read FFT size is smaller than acquired data.'); end
    if target(2) < sz(2), error('buildimages:phaseFftSizeTooSmall', 'Phase FFT size is smaller than acquired data.'); end
    if ~isfinite(target(3)) || target(3) < sz(3)
        error('buildimages:partitionFftSizeTooSmall', 'Partition FFT size is smaller than acquired data.');
    end
    mode = "ACS/Walsh";
    if ~isempty(scan.combination), mode = string(scan.combination); end
    bias = strcmpi(scan.biascorrect, 'ACS') && nSel > 1;
    needMaps = nSel > 1 && (mode == "ACS/Walsh" || bias);
    maps = [];
    if needMaps
        % Preserve mean(mean(...,field,'omitnan'),time,'omitnan').
        refSize = [sz(1:3) 1 1 nSel];
        total = zeros(refSize,'like',source);
        timeCount = zeros(refSize,'single');
        for t = 1:sz(4)
            fieldSum = zeros(refSize,'like',source);
            fieldCount = zeros(refSize,'single');
            for f = 1:sz(5)
                block = ffc_apply_whitening(source(:,:,:,t,f,:), info, sel, nc);
                valid = ~isnan(block); block(~valid) = 0;
                fieldSum = fieldSum + block;
                fieldCount = fieldCount + single(valid);
            end
            fieldMean = fieldSum ./ fieldCount;
            valid = ~isnan(fieldMean); fieldMean(~valid) = 0;
            total = total + fieldMean; timeCount = timeCount + single(valid);
        end
        reference = total ./ timeCount;
        clear total timeCount fieldSum fieldCount fieldMean block valid
        reference = windowkspace(reference, scan.window_size, scan.window_function);
        reference = ffc_center_pad(reference, target);
        maps = estimate_maps_from_acs_walsh(reference, 1:nSel, scan);
        clear reference
    end

    outSize = [target sz(4:5)];
    scan.magimage = zeros(outSize, class(source));
    scan.complexcombined = zeros(outSize, 'like', source);
    scan.magkspace = zeros(sz(1:5), class(source));
    keepCoils = nSel > 1 && mode ~= "Adaptive" && mode ~= "ACS/Walsh" && ...
        scan.retain_coil_images;
    scan.compleximage = [];
    if keepCoils
        scan.compleximage = zeros([outSize nSel], 'like', source);
    end
    ccopts = struct('coils',1:nSel,'sigma',4,'maskFrac',0.1,'referenceMode','signalcoil');
    if whiten, ccopts.referenceMode = 'rss'; end
    for f = 1:sz(5)
        for t = 1:sz(4)
            block = ffc_apply_whitening(source(:,:,:,t,f,:), info, sel, nc);
            scan.magkspace(:,:,:,t,f) = rssq(block,6);
            block = windowkspace(block, scan.window_size, scan.window_function);
            block = ffc_center_pad(block, target);
            if scan.TwoDimensional == 1
                coils = ifft2c(block);
            else
                coils = ifft3c(block);
            end
            clear block
            if nSel == 1
                combined = coils; magnitude = abs(combined);
            elseif mode == "Adaptive"
                combined = combine_adaptive_smoothed(coils,ccopts);
                magnitude = abs(combined);
            elseif mode == "ACS/Walsh"
                combined = combine_with_maps(coils,maps,1:nSel,scan);
                magnitude = abs(combined);
            else
                combined = sum(coils,6); magnitude = rssq(coils,6);
            end
            scan.complexcombined(:,:,:,t,f) = combined;
            scan.magimage(:,:,:,t,f) = magnitude;
            if keepCoils, scan.compleximage(:,:,:,t,f,:) = coils; end
        end
    end
    if ~keepCoils, scan.compleximage = scan.complexcombined; end
    if bias, scan.magimage = map_guided_bias_correct(scan.magimage,maps,scan); end
    scan.phaseimage = angle(scan.complexcombined);
end
