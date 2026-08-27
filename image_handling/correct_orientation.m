function obj = correct_orientation(obj)
%CORRECT_ORIENTATION Apply the prescribed FOV-offset phase ramps.
%
% The first three k-space dimensions are read, phase and partition. For a
% multislice 2D acquisition, dimension 3 indexes independent slices rather
% than kz and must not receive a through-plane Fourier phase ramp.

    if isempty(obj.complexkspace)
        return;
    end

    off1 = get_numeric_parameter(obj.param, 'OFF_CENTER_FIELD_OF_VIEW_1D', 0);
    off2 = get_numeric_parameter(obj.param, 'OFF_CENTER_FIELD_OF_VIEW_2D', 0);
    off3 = get_numeric_parameter(obj.param, 'OFF_CENTER_FIELD_OF_VIEW_3D', 0);

    fov1 = get_numeric_parameter(obj.param, 'FIELD_OF_VIEW', NaN);
    fov2 = get_numeric_parameter(obj.param, 'FIELD_OF_VIEW_PHASE', fov1);
    fov3 = get_numeric_parameter(obj.param, 'FIELD_OF_VIEW_3D', NaN);

    % Older formats do not always carry the Cameleon FOV fields. Their
    % ImageReconCore properties are stored in mm, whereas Cameleon offsets
    % and FOV values are stored in metres.
    if ~is_positive_finite(fov1) && ~isempty(obj.fov)
        fov1 = double(obj.fov) / 1000;
    end
    if ~is_positive_finite(fov2)
        if ~isempty(obj.fov_phase)
            fov2 = double(obj.fov_phase) / 1000;
        else
            fov2 = fov1;
        end
    end
    if ~is_positive_finite(fov3) && ~isempty(obj.fov_3d)
        fov3 = double(obj.fov_3d) / 1000;
    end

    nRead = size(obj.complexkspace, 1);
    nPhase = size(obj.complexkspace, 2);
    nPart = size(obj.complexkspace, 3);

    kRead = centred_frequency_axis(nRead, fov1, off1);
    kPhase = centred_frequency_axis(nPhase, fov2, off2);

    is3D = ~isempty(obj.TwoDimensional) && obj.TwoDimensional ~= 1;
    if is3D
        kPart = centred_frequency_axis(nPart, fov3, off3);
    else
        % In 2D multislice data this dimension is a slice label, not kz.
        kPart = zeros(1, nPart);
        off3 = 0;
    end

    [KRead, KPhase, KPart] = ndgrid(kRead, kPhase, kPart);
    phaseArgument = KRead .* off1 + KPhase .* off2 + KPart .* off3;
    phaseRamp = exp(-1i * 2*pi * phaseArgument);

    % Add singleton trailing dimensions so implicit expansion never touches
    % time, field or receiver axes.
    phaseRamp = reshape(phaseRamp, [nRead, nPhase, nPart, 1, 1, 1]);
    obj.complexkspace = obj.complexkspace .* cast(phaseRamp, ...
        'like', obj.complexkspace);
end


function k = centred_frequency_axis(n, fov, offset)
% FFT-shifted frequency coordinates, correct for both odd and even sizes.

    if offset == 0 || ~is_positive_finite(fov)
        k = zeros(1, n);
        return;
    end

    k = (-floor(n/2):ceil(n/2)-1) / fov;
end


function tf = is_positive_finite(value)
    tf = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end


function value = get_numeric_parameter(parameters, name, defaultValue)
% Read scalar numeric parameters from Cameleon structs or legacy cell lists.

    value = defaultValue;

    if isstruct(parameters) && isfield(parameters, name)
        candidate = first_numeric_scalar(parameters.(name));
        if ~isempty(candidate)
            value = candidate;
        end
        return;
    end

    if ~iscell(parameters) || isempty(parameters)
        return;
    end

    for row = 1:size(parameters, 1)
        key = parameters{row, 1};
        try
            keyMatches = strcmpi(strtrim(char(string(key))), name);
        catch
            keyMatches = false;
        end
        if ~keyMatches
            continue;
        end

        for col = 2:size(parameters, 2)
            candidate = first_numeric_scalar(parameters{row, col});
            if ~isempty(candidate)
                value = candidate;
                return;
            end
        end
    end
end


function value = first_numeric_scalar(candidate)
    value = [];

    if isnumeric(candidate) || islogical(candidate)
        candidate = double(candidate(:));
        candidate = candidate(isfinite(candidate));
        if ~isempty(candidate)
            value = candidate(1);
        end
        return;
    end

    if iscell(candidate)
        for index = 1:numel(candidate)
            value = first_numeric_scalar(candidate{index});
            if ~isempty(value)
                return;
            end
        end
        return;
    end

    try
        parsed = str2double(string(candidate));
        parsed = parsed(isfinite(parsed));
        if ~isempty(parsed)
            value = double(parsed(1));
        end
    catch
        value = [];
    end
end
