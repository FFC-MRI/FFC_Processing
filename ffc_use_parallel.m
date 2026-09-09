function tf = ffc_use_parallel(options, nJobs)
%FFC_USE_PARALLEL Never start a process pool as a side effect of interaction.
% Call parpool explicitly for large foreground batches. Workers run serially.
    tf = false;
    if nJobs < 16 || (isstruct(options) && isfield(options,'use_parallel') && ...
            ~options.use_parallel), return; end
    if ~license('test','Distrib_Computing_Toolbox'), return; end
    try
        if ~isempty(getCurrentTask()), return; end
        tf = ~isempty(gcp('nocreate'));
    catch
        tf = false;
    end
end
