function app = ffc_start()
%FFC_START Launch the editable performance edition with CVI contour support.
    app = ffc_start_cvi();
    if nargout == 0, clear app; end
end
