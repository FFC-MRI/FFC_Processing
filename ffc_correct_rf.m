function [out, info] = ffc_correct_rf(source, config, prior)
%FFC_CORRECT_RF Correct only acquired Cartesian lines of immutable k-space.
out = source;
if ~strcmp(config.status,'ready')
    info = struct('status',config.status,'version',config.version);
    return
end
if nargin<3, prior = []; end
if ~isempty(prior) && ~strcmp(prior.status,'corrected')
    info = prior;
    return
end
first = config.phaseRange(1);
last = config.phaseRange(2);
phase = first:last;
fullPhase = first==1 && last==size(source,2);
if fullPhase
    measured = source;
else
    measured = source(:,phase,:,:,:,:);
end
[corrected, info] = rs2d_predict_coherent_rf(measured,config.phaseOrder,config.is2D,prior);
info.phaseRange = config.phaseRange;
if strcmp(info.status,'corrected')
    if fullPhase
        out = corrected;
    else
        out(:,phase,:,:,:,:) = corrected;
    end
end
end
