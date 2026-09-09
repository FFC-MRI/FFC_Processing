function app = ffc_start_cvi()
%FFC_START_CVI Launch the CVI-contour test build.
%
% This build deliberately uses the uniquely named FFCProcessingCVI class.
% It therefore remains safe even if an older FFCProcessing.mlapp is still
% present in the same installation folder from a previous extraction.

    projectRoot = fileparts(mfilename('fullpath'));
    oldFolder = pwd;
    restoreFolder = onCleanup(@() cd(oldFolder)); %#ok<NASGU>
    cd(projectRoot);
    ffc_setup();

    try
        clear FFCProcessingCVI
    catch
    end
    rehash path

    resolved = which('FFCProcessingCVI');
    expected = fullfile(projectRoot,'FFCProcessingCVI.m');
    if isempty(resolved) || ~strcmpi(char(resolved),char(expected))
        error('ffc_start_cvi:WrongFFCProcessingCVI', ...
            ['MATLAB resolved FFCProcessingCVI to:\n  %s\n\nExpected:\n  %s\n' ...
             'Make sure this FFC_Processing folder is the current build and run ffc_start_cvi again.'], ...
            resolved,expected);
    end

    fprintf('Launching CVI contour test build 9 from:\n  %s\n',resolved);
    app = FFCProcessingCVI();
    if nargout == 0
        clear app
    end
end
