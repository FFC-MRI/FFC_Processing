function config = ffc_rf_correction_config(scan)
%FFC_RF_CORRECTION_CONFIG Standard RF stage policy and cache signature.
enabled = true;
K = scan.originalcomplexkspace;
ny = size(K,2);
config = struct('version',3,'status','ready','phaseOrder','auto', ...
                'is2D',isequal(double(scan.TwoDimensional),1) || size(K,3)==1, ...
                'phaseRange',[1 ny]);
if ~enabled
    config.status = 'disabled';
elseif ~isempty(scan.is_svs) && scan.is_svs
    config.status = 'spectroscopy';
elseif isempty(K) || size(K,1)<32 || ny<48
    config.status = 'insufficient_imaging_samples';
end
if ~strcmp(config.status,'ready'), return; end

p = scan.param;
plugin = '';
seq = '';
if isstruct(p)
    if isfield(p,'TRANSFORM_PLUGIN'), plugin = char(p.TRANSFORM_PLUGIN); end
    if isfield(p,'SEQUENCE_NAME'), seq = char(p.SEQUENCE_NAME); end
end
if ~isempty(regexpi([plugin ' ' seq ' ' char(scan.sequence)], ...
        'spiral|radial|propeller|stackofstars|stack.of.stars|non.?cartesian|blade', 'once')) ...
        || ~isempty(scan.bladeangles) && numel(scan.bladeangles)>1
    config.status = 'unsupported_trajectory';
    return
end
if strncmpi(plugin,'Sequential',10) || strcmpi(plugin,'linear')
    config.phaseOrder = 'linear';
elseif strcmpi(plugin,'CentricZigZag2D')
    config.phaseOrder = 'centric-auto';
elseif strcmpi(plugin,'Centric2D')
    config.phaseOrder = 'centric2d';
end

% Preprocessing centre-pads acquired phase lines before reconstruction.
% Fit only that original span, then put the estimate back at its old indices.
nAcquired = ny;
if ~isempty(scan.rawdata)
    nAcquired = min(ny,size(scan.rawdata,2));
elseif isstruct(p) && isfield(p,'MATRIX_DIMENSION_2D')
    stored = p.MATRIX_DIMENSION_2D;
    if isnumeric(stored) && isscalar(stored) && isfinite(stored) && stored>=1 && fix(stored)==stored
        nAcquired = min(ny,double(stored));
    end
end
start = floor((ny-nAcquired)/2)+1;
config.phaseRange = [start start+nAcquired-1];
if nAcquired<48, config.status = 'insufficient_acquired_phase_lines'; end
end
