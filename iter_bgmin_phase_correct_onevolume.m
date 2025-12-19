function [kspaceOut, out] = iter_bgmin_phase_correct_onevolume(kspaceIn, dims, opts)
%ITER_BGMIN_PHASE_CORRECT_ONEVOLUME
% Correct ONLY one volume: slice=7, bevo=end, tevo=end (or as specified).
%
% dims: struct specifying which dimension numbers correspond to slice/bevo/tevo
%   dims.sliceDim = ...
%   dims.bevoDim  = ...
%   dims.tevoDim  = ...
%
% opts: passed to iter_bgmin_phase_correct (your iterative background minimiser)
%   opts.sliceIndex (default 7)
%   opts.bevoIndex  (default 'end')
%   opts.tevoIndex  (default 'end')

arguments
    kspaceIn {mustBeNumeric}
    dims struct
    opts struct = struct()
end

% defaults
if ~isfield(opts,'sliceIndex'), opts.sliceIndex = 7; end
if ~isfield(opts,'bevoIndex'),  opts.bevoIndex  = 'end'; end
if ~isfield(opts,'tevoIndex'),  opts.tevoIndex  = 'end'; end

sz = size(kspaceIn);

% Resolve 'end'
sliceIdx = resolveIndex(opts.sliceIndex, sz(dims.sliceDim));
bevoIdx  = resolveIndex(opts.bevoIndex,  sz(dims.bevoDim));
tevoIdx  = resolveIndex(opts.tevoIndex,  sz(dims.tevoDim));

% Build indexing cell
subs = repmat({':'}, 1, ndims(kspaceIn));
subs{dims.sliceDim} = sliceIdx;
subs{dims.bevoDim}  = bevoIdx;
subs{dims.tevoDim}  = tevoIdx;

% Extract just that subset (keeps kx, ky, coils, etc.)
kSub = kspaceIn(subs{:});

% Run your iterative correction on the subset only
% (This calls the optimizer you already have: iter_bgmin_phase_correct)
[kSubCorr, out] = correct_ky_parity_phase(kSub, opts);

% Insert back
kspaceOut = kspaceIn;
kspaceOut(subs{:}) = kSubCorr;

end

function idx = resolveIndex(val, dimLen)
if ischar(val) || (isstring(val) && val=="end")
    idx = dimLen;
else
    idx = val;
end
end
