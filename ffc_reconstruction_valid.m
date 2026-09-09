function valid = ffc_reconstruction_valid(scan)
%FFC_RECONSTRUCTION_VALID Display changes and postprocessing preserve validity.
    valid = ~isempty(scan.magimage) && ~isempty(scan.reconstruction_key) && ...
        isequaln(scan.reconstruction_key, ffc_reconstruction_key(scan));
end
