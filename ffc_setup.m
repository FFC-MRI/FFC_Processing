function projectRoot = ffc_setup()
%FFC_SETUP Add the active pipeline and retain legacy sequence deserialisation.
% Use this instead of adding every experimental/demo directory recursively.
    projectRoot = fileparts(mfilename('fullpath'));
    folders = {'misc_functions','file_handling','file_handling/Sequence files', ...
        'file_handling/wave generators','kspace_handling','image_handling', ...
        'analysis_methods','contour_handling','analysis_methods/t1_mapping', ...
        'analysis_methods/t1_mapping/mfiles','pipeline_methods'};
    for i = numel(folders):-1:1
        folder = fullfile(projectRoot,folders{i});
        if isfolder(folder), addpath(folder,'-begin'); end
    end
    % Canonical filter defaults take priority over the old duplicate wrapper.
    addpath(fullfile(projectRoot,'image_handling'),'-begin');
    addpath(projectRoot,'-begin');
end
