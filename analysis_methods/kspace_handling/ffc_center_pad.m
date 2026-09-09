function out = ffc_center_pad(in, target)
%FFC_CENTER_PAD Centre spatial zero filling with a single final allocation.
% Odd padding follows the established floor(pre), ceil(post) convention.
    sz = size(in); sz(end+1:3) = 1;
    target = double(target(:).');
    if numel(target) ~= 3 || any(~isfinite(target)) || ...
            any(target < sz(1:3)) || any(target ~= round(target))
        error('ffc_center_pad:InvalidSize', 'Target must contain three integer dimensions at least as large as the input.');
    end
    delta = target - sz(1:3);
    if all(delta == 0), out = in; return; end
    pre = floor(delta / 2);
    out = zeros([target sz(4:end)], 'like', in);
    ix = repmat({':'}, 1, max(3, ndims(in)));
    for d = 1:3, ix{d} = pre(d) + (1:sz(d)); end
    out(ix{:}) = in;
end
