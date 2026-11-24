function [Kout, info] = ffc_tone_prototype_hook(Kin, config)
%FFC_TONE_PROTOTYPE_HOOK Optional background-validated RF subtraction.
% Called after preprocessing and before any FOV phase, window or whitening.
% No filesystem writes or persistent state. Always pass immutable source data.
if nargin<2, config = ffc_tone_prototype_config(); end
if config.enabled
    [Kout,info] = rs2d_predict_coherent_rf(Kin,config.phaseOrder);
else
    Kout = Kin;
    info = struct('status','disabled_or_ineligible','version',config.version);
end
end
