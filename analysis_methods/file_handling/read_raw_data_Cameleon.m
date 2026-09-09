function [data, par, data_name] = read_raw_data_Cameleon(data_path, options)
% More efficient version:
%   - Parse XML header once into par.cameleon
%   - Avoid repeated DOM traversals and str2num calls

arguments
    data_path (1, :) char = uigetdir(pwd) % if not specified, use the default folder
    options = []
end

%----------------------
%-- Parameters
%----------------------
if data_path == 0  % uigetdir does that if you hit 'cancel'
    error("Data read canceled by user.");
end

[~, data_name, extension] = fileparts(data_path);
folder_names = regexp(data_path, filesep, 'split'); %#ok<NASGU>
data_name = [data_name extension];

headerFile = [data_path filesep 'header.xml'];
doc = xmlread(headerFile);

%--------------------------------------------------
% Parse ALL parameters once into par.cameleon
%--------------------------------------------------
par = struct();
par.cameleon = parseCameleonParameters(doc);

% Convenience numeric getter
getNum = @(name, default) getNumericFromStruct(par.cameleon, name, default);

% /!\ If not in this list, all the sequence parameters are always
% accessible as par.cameleon.<PARAMETER_NAME>

% Dimensions & acquisition settings
par.md1d             = getNum('MATRIX_DIMENSION_1D',              0);
par.md2d             = getNum('MATRIX_DIMENSION_2D',              0);
par.md3d             = getNum('MATRIX_DIMENSION_3D',              0);
par.md4d             = getNum('MATRIX_DIMENSION_4D',              0);
par.md5d             = getNum('USER_MATRIX_DIMENSION_5D',         0);
par.accudim          = getNum('ACCU_DIM',                         0);
par.ncoils           = getNum('RECEIVER_COUNT',                   1);
par.navg             = getNum('NUMBER_OF_AVERAGES',               1);
par.etl              = getNum('ECHO_TRAIN_LENGTH',                0);
par.sw               = getNum('SPECTRAL_WIDTH',                   0);
par.te               = getNum('ECHO_TIME',                        0);
par.ti_step          = getNum('TI_STEP',                          0); % for T1 IR
par.dte              = getNum('DELTA_TE',                         0); % for B0 mapping
par.tx_step          = getNum('TX_DELTA_STEPS',                   0); % for RF calibration
par.resx             = getNum('RES_X',                            0);
par.resy             = getNum('RES_Y',                            0);
par.resz             = getNum('RES_Z',                            0);
par.amd1d            = getNum('ACQUISITION_MATRIX_DIMENSION_1D',  0);
par.amd2d            = getNum('ACQUISITION_MATRIX_DIMENSION_2D',  0);
par.amd3d            = getNum('ACQUISITION_MATRIX_DIMENSION_3D',  0);
par.amd4d            = getNum('ACQUISITION_MATRIX_DIMENSION_4D',  0);
par.amd5d            = getNum('ACQUISITION_MATRIX_DIMENSION_5D',  0);
par.umd1d            = getNum('USER_MATRIX_DIMENSION_1D',         0);
par.umd2d            = getNum('USER_MATRIX_DIMENSION_2D',         0);
par.umd3d            = getNum('USER_MATRIX_DIMENSION_3D',         0);
par.umd4d            = getNum('USER_MATRIX_DIMENSION_4D',         0);
par.nbgsteps         = getNum('NUMBER_GRADIENT_STEPS',            0); % gradient calib
par.gampmin          = getNum('GRADIENT_AMP_MIN',                 0);
par.gampmax          = getNum('GRADIENT_AMP_MAX',                 0);
par.partFourierDim2D = getNum('USER_PARTIAL_FOURIER_2D',          0); % partial k-space
par.seqTime          = getNum('SEQUENCE_TIME',                    0);
par.ringTime         = getNum('TX_TIME_RINGING',                  0);
par.tiFirstDelay     = getNum('TI_FIRST_DELAY',                   0);

par.fa               = getNum('FLIP_ANGLE',                       0);
par.tx90length       = getNum('TX_90PULSE_LENGTH',                0);
par.tx_shape_max     = getNum('TX_SHAPE_MAX',                     0);

% Sequence name/version (string)
if isfield(par.cameleon, 'SEQUENCE_NAME')
    par.seqName = char(par.cameleon.SEQUENCE_NAME);
else
    par.seqName = '';
end
par.serieName = getSerieNameFromSerieXML(data_path);
if isfield(par.cameleon, 'SEQUENCE_VERSION')
    par.seqVersion = char(par.cameleon.SEQUENCE_VERSION);
else
    par.seqVersion = '';
end

%----------------------
%-- Read Data
%----------------------
dataFile = [data_path filesep 'data.dat'];

% 'ieee-be' = big endian (file written by Java)
nComplex = par.md1d * par.md2d * par.md3d * par.md4d * par.ncoils;
fid = fopen(dataFile, 'rb');
if fid < 0
    error('Could not open data file: %s', dataFile);
end

data = fread(fid, 2 * nComplex, 'float32=>single', 0, 'ieee-be');
fclose(fid);

% Reshape and form complex data: convention [real; imag] or [imag; real]?
% Original code: data = complex(data(2,:,:,:,:,:), data(1,:,:,:,:,:));
data = reshape(data, [2, par.md1d, par.md2d, par.md3d, par.md4d, par.ncoils]);
data = complex(data(2,:,:,:,:,:), data(1,:,:,:,:,:));

z = size(data);
if numel(z) > 2
    data = reshape(data, z(2:end));
else
    % 1D data case
    if (par.etl ~= 0)
        data = reshape(data, z(2)/par.etl, par.etl);
    else
        data = data.'; % transpose so points are along first dim
    end
end

%------------------------------------------------------
% Optional sequence-specific processing (MRF)
%------------------------------------------------------
if contains(par.seqName, 'MRF')  % process data from a MRF sequence
    if par.accudim == 4
        % MRF v2
        test_mrf_v2 = 1;
        test_mrf_v3 = 0;
    else
        test_mrf_v2 = 0;
        if size(data,4)/par.navg > 1
            test_mrf_v3 = 1;
        else
            test_mrf_v3 = 0;
        end
    end

    if test_mrf_v3 == 1
        % MRF v3: average across navg, more efficient preallocation
        sz = size(data); % [x y z t]
        nTime = sz(4) / par.navg;

        if rem(sz(4), par.navg) ~= 0
            error('MRF v3: size(data,4) must be divisible by par.navg.');
        end

        % Reorder into [x y z nTime navg] by explicit indexing (still vectorised)
        data_r = zeros(sz(1), sz(2), sz(3), nTime, par.navg, 'like', data);
        for n = 1:par.navg
            data_r(:,:,:,:,n) = data(:,:,:,n:par.navg:end);
        end

        data = sum(data_r, 5);           % sum over navg
        data = permute(data, [1 3 4 2]); % same as original

    elseif test_mrf_v2 == 1
        % MRF v2: make RX channel dimension behave
        data = permute(data, [1 3 4 2 5]);
    end
end

end

%======================================================================
% Helper: parse all <entry> tags efficiently into a struct
%======================================================================
function cameleon = parseCameleonParameters(doc)
% Reads all key-value pairs from XML into a struct:
%   cameleon.<key> = numeric scalar/vector, boolean, or char

items = doc.getElementsByTagName('entry');
nItems = items.getLength;

cameleon = struct();

for k = 0:nItems-1
    item = items.item(k);

    % Key
    keyNode = item.getElementsByTagName('key').item(0);
    if isempty(keyNode) || isempty(keyNode.getFirstChild())
        continue;
    end
    parname = char(keyNode.getFirstChild.getData());

    % Values
    valueNode = item.getElementsByTagName('value').item(0);
    if isempty(valueNode)
        cameleon.(parname) = '';
        continue;
    end

    valNodes = valueNode.getElementsByTagName('value');
    nVals    = valNodes.getLength;

    if nVals == 0
        cameleon.(parname) = '';
        continue;
    end

    vals = cell(1, nVals);
    isBool = false;

    for j = 0:nVals-1
        vNode = valNodes.item(j);
        txtNode = vNode.getFirstChild();
        if isempty(txtNode)
            vals{j+1} = '';
            continue;
        end
        vStr = char(txtNode.getData());

        % Boolean?
        if strcmp(vStr, 'true')
            vals{j+1} = true;
            isBool = true;
        elseif strcmp(vStr, 'false')
            vals{j+1} = false;
            isBool = true;
        else
            vals{j+1} = strtrim(vStr);
        end
    end

    % Reduce cell array to appropriate type
    if isBool
        % If any boolean values, store as logical array or scalar
        if nVals == 1
            cameleon.(parname) = logical(vals{1});
        else
            cameleon.(parname) = logical([vals{:}]);
        end
    else
        % Concatenate as a single string
        flat = strtrim(strjoin(vals(~cellfun(@isempty, vals)), ' '));

        if isempty(flat)
            cameleon.(parname) = '';
        else
            % Try fast scalar numeric parse first
            numScalar = str2double(flat);
            if ~isnan(numScalar) && ~contains(flat, ' ')
                cameleon.(parname) = numScalar;
            else
                % As a fallback, try vector parse once (still much cheaper than
                % calling str2num many times in loops)
                numVec = str2num(flat); %#ok<ST2NM>
                if ~isempty(numVec) && isnumeric(numVec)
                    cameleon.(parname) = numVec;
                else
                    cameleon.(parname) = flat;
                end
            end
        end
    end
end
end

%======================================================================
% Helper: safe numeric getter from a struct
%======================================================================
function val = getNumericFromStruct(s, fieldName, default)
% Returns numeric scalar from s.(fieldName) if possible, otherwise default.

if nargin < 3
    default = 0;
end

if ~isfield(s, fieldName)
    val = default;
    return;
end

raw = s.(fieldName);

if isnumeric(raw)
    if isempty(raw)
        val = default;
    else
        val = raw(1); % in case it's a vector, take first element
    end
elseif islogical(raw)
    val = double(raw(1));
else
    % String or char: try to parse
    numVal = str2double(string(raw));
    if isnan(numVal)
        val = default;
    else
        val = numVal;
    end
end

end
function serieName = getSerieNameFromSerieXML(data_path)
%GETSERIENAMFROMSERIEXML Extracts <name>...</name> from Serie.xml in data_path.
serieName = '';

seriePath = fullfile(data_path, 'Serie.xml');
if ~isfile(seriePath)
    return;
end

try
    txt = fileread(seriePath);

    % First <name> tag inside <serie>
    tok = regexp(txt, '<name>\s*(.*?)\s*</name>', 'tokens', 'once');
    if ~isempty(tok)
        serieName = strtrim(tok{1});

        % Minimal XML entity decoding
        serieName = strrep(serieName, '&amp;', '&');
        serieName = strrep(serieName, '&lt;',  '<');
        serieName = strrep(serieName, '&gt;',  '>');
        serieName = strrep(serieName, '&quot;','"');
        serieName = strrep(serieName, '&apos;','''');
    end
catch
    % leave as ''
end
end
