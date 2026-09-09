function scan = ffc_invalidate_analysis(scan)
%FFC_INVALIDATE_ANALYSIS Derived fits must follow the current magnitude data.
    scan.T1Maps = [];
    scan.R1Maps = [];
    scan.T1FitResults = [];
    scan.dispersioncurve = [];
    if isstruct(scan.param)
        fields = fieldnames(scan.param);
        scan.param = rmfield(scan.param,fields(startsWith(fields,'SimpleIR')));
    end
    scan.magnitude_generation = scan.magnitude_generation + 1;
end
