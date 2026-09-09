function geometry = cvi_dicom_frame_geometry(info,frameNumber)
%CVI_DICOM_FRAME_GEOMETRY Geometry and frame indices of one DICOM frame.
%
% In addition to the image plane geometry, enhanced-MR FrameContentSequence
% metadata are returned when available.  This avoids assuming a particular
% slice/phase frame ordering in the multi-frame object.

    if nargin < 2, frameNumber = 1; end
    frame = struct();
    if isfield(info,'PerFrameFunctionalGroupsSequence')
        seq = info.PerFrameFunctionalGroupsSequence;
        itemName = sprintf('Item_%d',frameNumber);
        if isfield(seq,itemName)
            frame = seq.(itemName);
        else
            names = fieldnames(seq);
            if frameNumber <= numel(names), frame = seq.(names{frameNumber}); end
        end
    end
    shared = struct();
    if isfield(info,'SharedFunctionalGroupsSequence')
        shared = firstSequenceItem(info.SharedFunctionalGroupsSequence);
    end

    ipp = firstAvailableValue(frame,shared,info, ...
        {'PlanePositionSequence','PlanePositionPatientSequence'},'ImagePositionPatient');
    iop = firstAvailableValue(frame,shared,info, ...
        {'PlaneOrientationSequence','PlaneOrientationPatientSequence'},'ImageOrientationPatient');
    spacing = firstAvailableValue(frame,shared,info, ...
        {'PixelMeasuresSequence'},'PixelSpacing');

    ipp = double(ipp(:)); iop = double(iop(:)); spacing = double(spacing(:));
    if numel(ipp) ~= 3 || numel(iop) ~= 6 || numel(spacing) < 2 || ...
            any(~isfinite([ipp;iop;spacing(1:2)]))
        error('cvi_dicom_frame_geometry:MissingGeometry', ...
            'Could not recover ImagePositionPatient/ImageOrientationPatient/PixelSpacing for frame %d.',frameNumber);
    end
    row = iop(1:3); row = row ./ norm(row);
    column = iop(4:6); column = column - row .* dot(row,column);
    column = column ./ norm(column);

    % Enhanced DICOM frame ordering.  Philips cardiac cine in the supplied
    % study has DimensionIndexValues = [stack, inStackSlice, cardiacPhase].
    frameContent = struct();
    if isstruct(frame) && isfield(frame,'FrameContentSequence')
        frameContent = firstSequenceItem(frame.FrameContentSequence);
    end
    inStack = scalarField(frameContent,'InStackPositionNumber',NaN);
    temporal = scalarField(frameContent,'TemporalPositionIndex',NaN);
    dimValues = vectorField(frameContent,'DimensionIndexValues',[]);

    phaseIndex = NaN;
    if ~isempty(dimValues)
        % Last dimension is the varying cardiac-phase dimension for the
        % Philips enhanced cine used here.  Validate later against the inferred
        % number of phases before accepting it.
        phaseIndex = dimValues(end);
    elseif isfinite(temporal) && temporal > 1
        phaseIndex = temporal;
    end

    geometry = struct();
    geometry.imagePositionPatient = ipp;
    geometry.rowDirection = row;
    geometry.columnDirection = column;
    geometry.normalDirection = cross(row,column);
    geometry.rowSpacingMm = spacing(1);
    geometry.columnSpacingMm = spacing(2);
    geometry.inStackPositionNumber = inStack;
    geometry.temporalPositionIndex = temporal;
    geometry.dimensionIndexValues = dimValues;
    geometry.phaseIndex = phaseIndex;
end

function value = firstAvailableValue(frame,shared,info,sequences,field)
    value = [];
    parents = {frame,shared};
    for p = 1:numel(parents)
        parent = parents{p};
        for s = 1:numel(sequences)
            if isstruct(parent) && isfield(parent,sequences{s})
                item = firstSequenceItem(parent.(sequences{s}));
                if isfield(item,field) && ~isempty(item.(field))
                    value = item.(field); return;
                end
            end
        end
        if isstruct(parent) && isfield(parent,field) && ~isempty(parent.(field))
            value = parent.(field); return;
        end
    end
    if isfield(info,field) && ~isempty(info.(field)), value = info.(field); end
end

function item = firstSequenceItem(sequence)
    item = struct();
    if ~isstruct(sequence), return; end
    names = fieldnames(sequence);
    if ~isempty(names), item = sequence.(names{1}); end
end

function value = scalarField(s,name,defaultValue)
    value = defaultValue;
    if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name)), return; end
    candidate = double(s.(name));
    candidate = candidate(:);
    candidate = candidate(isfinite(candidate));
    if ~isempty(candidate), value = candidate(1); end
end

function value = vectorField(s,name,defaultValue)
    value = defaultValue;
    if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name)), return; end
    candidate = double(s.(name));
    candidate = candidate(:).';
    candidate = candidate(isfinite(candidate));
    if ~isempty(candidate), value = candidate; end
end
