function mag_corr = apply_bias_field(mag, bias, mask3, obj)
% Applies mag./bias only inside mask; outside unchanged.
% Also clamps bias to avoid extreme amplification.

    sz = size(mag); sz(end+1:5)=1;
    X=sz(1); Y=sz(2); S=sz(3); T=sz(4); F=sz(5);
    mag = reshape(mag, [X Y S T F]);

    % Clamp settings (good defaults)
    bmin = 0.4;  % don't allow bias < 0.4 (would amplify >2.5x)
    bmax = 2.5;  % don't allow bias > 2.5
    if isfield(obj,'bias_clip_min') && ~isempty(obj.bias_clip_min), bmin = obj.bias_clip_min; end
    if isfield(obj,'bias_clip_max') && ~isempty(obj.bias_clip_max), bmax = obj.bias_clip_max; end

    bias = single(bias);
    bias = max(bmin, min(bmax, bias));        % clamp
    mask3 = logical(mask3);

    mag_corr = mag; % default unchanged

    for s = 1:S
        m = mask3(:,:,s);
        if ~any(m(:)), continue; end

        b = bias(:,:,s);
        for t = 1:T
            for f = 1:F
                I = mag(:,:,s,t,f);
                Icorr = I;
                Icorr(m) = I(m) ./ (b(m) + eps('single'));
                mag_corr(:,:,s,t,f) = Icorr;
            end
        end
    end
end
