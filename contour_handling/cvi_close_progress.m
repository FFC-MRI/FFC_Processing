function cvi_close_progress(progress)
%CVI_CLOSE_PROGRESS Safely close a progress dialog used by CVI import.
    try
        if ~isempty(progress) && isvalid(progress)
            close(progress);
        end
    catch
        % Import error handling must not be masked by dialog cleanup.
    end
end
