function scan = ffc_normalize_viewer_scan(scan)
%FFC_NORMALIZE_VIEWER_SCAN Apply viewer defaults BEFORE validating/rebuilding.
% Idempotent. Numerical core callers keep their per-coil output policy.
    scan.retain_coil_images = false;
    scan.retain_phase_images = false;
    if isempty(scan.fft_size), scan.fft_size = scan.samples; end
    if isempty(scan.fft_size_3d), scan.fft_size_3d = scan.slices; end
    if isempty(scan.window_size), scan.window_size = 1; end
    choices = {'None','Hann','Blackman','Hamming','Kaiser-Bessel','Tukey'};
    match = find(strcmpi(scan.window_function,choices),1);
    if isempty(match), match = 1; end
    scan.window_function = choices{match};
    n = max(1,double(scan.n_receivers));
    sel = scan.multichannel_recon;
    if islogical(sel), sel = find(sel(:).'); end
    sel = unique(round(double(sel(:).')), 'stable');
    sel = sel(isfinite(sel) & sel >= 1 & sel <= n);
    if isempty(sel), sel = 1:n; end
    scan.multichannel_recon = sel;
    if isempty(scan.combination) || ~any(strcmp(scan.combination,{'RSSQ','Adaptive','ACS/Walsh'}))
        scan.combination = 'Adaptive';
    end
    if isempty(scan.biascorrect) || ~any(strcmp(scan.biascorrect,{'Off','ACS'}))
        scan.biascorrect = 'Off';
    end
    if isempty(scan.noise_whitening), scan.noise_whitening = n > 1; end
    scan.noise_whitening = n > 1 && logical(scan.noise_whitening);
    if isempty(scan.denoise_filter), scan.denoise_filter = 'None'; end
    if isempty(scan.denoise_params) || ~isstruct(scan.denoise_params)
        scan.denoise_params = ffc_default_denoise_params(scan.denoise_filter);
    end
end
