function meta = ffc_extract_dicom_metadata(scan)
%FFC_EXTRACT_DICOM_METADATA Extract DICOM identity/study fields from FCI raw metadata.
%
% No DICOM source file is required.  The preferred source is the top-level
% metadata stored in Serie.xml (read into scan.source_metadata), followed by
% the complete Cameleon parameter structure in scan.param.

meta = struct( ...
    'PatientName','', ...
    'PatientID','', ...
    'PatientBirthDate','', ...
    'PatientSex','', ...
    'PatientWeight','', ...
    'StudyDescription','Field-Cycling MRI', ...
    'StudyID','', ...
    'StudyDate','', ...
    'StudyTime','', ...
    'AcquisitionDateTime','', ...
    'PatientPosition','', ...
    'SeriesDescription','', ...
    'SeriesID','', ...
    'AcquisitionUUID','', ...
    'ParticipantKey','', ...
    'StudyKey','', ...
    'MetadataSource','raw dataset');

sources = {};
% Re-read the raw XML context at export time as well as using the metadata
% retained at load time.  This catches patient/subject fields stored in
% participant-level XML or nested blocks that older loaders did not flatten.
rawXml = readRawXmlIdentity(scan);
if isstruct(rawXml) && ~isempty(fieldnames(rawXml))
    sources{end+1} = rawXml; %#ok<AGROW>
end
try
    if isobject(scan) && isprop(scan,'source_metadata') && isstruct(scan.source_metadata)
        sources{end+1} = scan.source_metadata; %#ok<AGROW>
    elseif isstruct(scan) && isfield(scan,'source_metadata') && isstruct(scan.source_metadata)
        sources{end+1} = scan.source_metadata; %#ok<AGROW>
    end
catch
end
try
    if isobject(scan) && isprop(scan,'param') && isstruct(scan.param)
        sources{end+1} = scan.param; %#ok<AGROW>
    elseif isstruct(scan) && isfield(scan,'param') && isstruct(scan.param)
        sources{end+1} = scan.param; %#ok<AGROW>
    end
catch
end

rawPatient = firstText(sources, { ...
    'PatientName','PATIENT_NAME','Patient','PATIENT','SubjectName','SUBJECT_NAME'});
meta.PatientName = formatPatientName(rawPatient);
meta.PatientID = firstText(sources, { ...
    'PatientID','PATIENT_ID','PatientIdentifier','PATIENT_IDENTIFIER', ...
    'SubjectID','SUBJECT_ID','SubjectNumber','SUBJECT_NUMBER','PatientNumber','PATIENT_NUMBER'});

meta.PatientBirthDate = normalizeDicomDate(firstText(sources, { ...
    'PatientBirthDate','PATIENT_BIRTH_DATE','PatientBirthday','PATIENT_BIRTHDAY', ...
    'PatientBirthDay','BirthDate','BIRTH_DATE','DateOfBirth','DATE_OF_BIRTH', ...
    'DOB','SubjectBirthDate','SUBJECT_BIRTH_DATE'}));
meta.PatientSex = normalizeSex(firstText(sources, { ...
    'PatientSex','PATIENT_SEX','Sex','SEX','Gender','GENDER', ...
    'SubjectSex','SUBJECT_SEX','SubjectGender','SUBJECT_GENDER'}));
meta.PatientWeight = firstText(sources, { ...
    'PatientWeight','PATIENT_WEIGHT','SubjectWeight','SUBJECT_WEIGHT','WeightKg','WEIGHT_KG'});

study = firstText(sources, {'StudyDescription','STUDY_DESCRIPTION','Study','STUDY','StudyName','STUDY_NAME'});
if ~isempty(study), meta.StudyDescription = study; end
meta.StudyID = firstText(sources, {'StudyID','STUDY_ID','StudyIdentifier','STUDY_IDENTIFIER','StudyNumber','STUDY_NUMBER'});
if isempty(meta.StudyID), meta.StudyID = meta.StudyDescription; end
session = firstText(sources, {'Session','SESSION','SessionDate','SESSION_DATE'});
meta.AcquisitionDateTime = firstText(sources, { ...
    'AcquisitionDateTime','ACQUISITION_DATETIME','AcquisitionDate','ACQUISITION_DATE'});
[acqDate,acqTime] = parseDateTime(meta.AcquisitionDateTime);
sessionDate = normalizeDicomDate(session);
if ~isempty(sessionDate), meta.StudyDate = sessionDate; else, meta.StudyDate = acqDate; end
meta.StudyTime = acqTime;

position = firstText(sources, {'PatientPosition','PATIENT_POSITION','SubjectPosition','SUBJECT_POSITION'});
meta.PatientPosition = normalizePatientPosition(position);
meta.SeriesDescription = firstText(sources, {'SerieName','SeriesName','SERIES_NAME','ProtocolName','PROTOCOL_NAME'});
if isempty(meta.SeriesDescription)
    meta.SeriesDescription = scanText(scan,'displayname');
end
if isempty(meta.SeriesDescription)
    meta.SeriesDescription = scanText(scan,'sequence');
end
meta.SeriesID = firstText(sources, {'SerieID','SeriesID','SERIES_ID'});
meta.AcquisitionUUID = firstText(sources, {'AcquisitionUUID','ACQUISITION_UUID','SerieUUID','SeriesUUID'});

% Stable matching keys used only internally.  They are never written as PHI.
participantBase = rawPatient;
if isempty(participantBase), participantBase = meta.PatientID; end
if isempty(participantBase), participantBase = meta.PatientName; end
meta.ParticipantKey = normalizeKey(participantBase);
studyBase = strjoin(nonEmpty({meta.ParticipantKey,meta.StudyDescription,meta.StudyDate}), '|');
meta.StudyKey = normalizeKey(studyBase);
end

function value = scanText(scan,name)
value = '';
try
    if (isobject(scan) && isprop(scan,name)) || (isstruct(scan) && isfield(scan,name))
        value = textValue(scan.(name));
    end
catch
end
end

function value = firstText(sources,names)
value = '';
for s = 1:numel(sources)
    src = sources{s};
    if ~isstruct(src), continue; end
    fields = fieldnames(src);
    normFields = cellfun(@normalizeFieldName,fields,'UniformOutput',false);
    for n = 1:numel(names)
        wanted = normalizeFieldName(names{n});
        idx = find(strcmp(normFields,wanted),1,'first');
        if isempty(idx), continue; end
        candidate = textValue(src.(fields{idx}));
        if ~isempty(candidate)
            value = candidate;
            return
        end
    end
end
end

function out = normalizeFieldName(in)
out = upper(regexprep(char(string(in)),'[^A-Za-z0-9]',''));
end

function out = textValue(value)
out = '';
if isempty(value), return; end
if iscell(value)
    if isempty(value), return; end
    value = value{1};
end
if isstring(value)
    value = char(value(1));
elseif isnumeric(value) || islogical(value)
    if isscalar(value) && isfinite(double(value))
        value = num2str(value);
    else
        return
    end
end
if ischar(value)
    out = strtrim(value);
end
end

function out = formatPatientName(raw)
out = strtrim(raw);
if isempty(out), return; end
if contains(out,'^'), return; end
parts = regexp(out,'\s+','split');
parts = parts(~cellfun(@isempty,parts));
if numel(parts) == 2
    out = [parts{1} '^' parts{2}];
elseif numel(parts) > 2
    out = [parts{1} '^' strjoin(parts(2:end),' ')];
end
end

function out = participantLabel(raw)
out = '';
raw = strtrim(raw);
if isempty(raw), return; end
parts = regexp(raw,'[\s\^]+','split');
parts = parts(~cellfun(@isempty,parts));
if ~isempty(parts), out = parts{1}; end
end

function out = normalizeSex(raw)
out = '';
r = lower(strtrim(raw));
if isempty(r), return; end
if any(strcmp(r,{'m','male','man'})), out = 'M'; return; end
if any(strcmp(r,{'f','female','woman'})), out = 'F'; return; end
if any(strcmp(r,{'o','other','unknown','u'})), out = 'O'; return; end
end

function out = normalizeDicomDate(raw)
out = '';
r = strtrim(raw);
if isempty(r), return; end
% YYYY-MM-DD / YYYYMMDD / ISO-8601
m = regexp(r,'(?<!\d)(\d{4})[-/.]?(\d{2})[-/.]?(\d{2})(?!\d)','tokens','once');
if ~isempty(m)
    out = [m{1} m{2} m{3}];
    return
end
% DD-MM-YYYY / DD/MM/YYYY
m = regexp(r,'(?<!\d)(\d{1,2})[-/.](\d{1,2})[-/.](\d{4})(?!\d)','tokens','once');
if ~isempty(m)
    out = sprintf('%04d%02d%02d',str2double(m{3}),str2double(m{2}),str2double(m{1}));
end
end

function [dateText,timeText] = parseDateTime(raw)
dateText = normalizeDicomDate(raw);
timeText = '';
r = strtrim(raw);
if isempty(r), return; end
m = regexp(r,'[T ](\d{2}):?(\d{2}):?(\d{2})(?:\.(\d+))?','tokens','once');
if isempty(m), return; end
timeText = [m{1} m{2} m{3}];
if numel(m) >= 4 && ~isempty(m{4})
    timeText = [timeText '.' m{4}];
end
end

function out = normalizePatientPosition(raw)
out = '';
r = upper(regexprep(strtrim(raw),'[^A-Za-z]',''));
switch r
    case {'HEADFIRSTSUPINE','HFS'}, out = 'HFS';
    case {'HEADFIRSTPRONE','HFP'}, out = 'HFP';
    case {'FEETFIRSTSUPINE','FFS'}, out = 'FFS';
    case {'FEETFIRSTPRONE','FFP'}, out = 'FFP';
    case {'HEADFIRSTDECUBITUSRIGHT','HFDR'}, out = 'HFDR';
    case {'HEADFIRSTDECUBITUSLEFT','HFDL'}, out = 'HFDL';
    case {'FEETFIRSTDECUBITUSRIGHT','FFDR'}, out = 'FFDR';
    case {'FEETFIRSTDECUBITUSLEFT','FFDL'}, out = 'FFDL';
end
end

function metadata = readRawXmlIdentity(scan)
%READRAWXMLIDENTITY Recover participant/study identity directly from raw XML.
% This intentionally has no DICOM dependency.  It searches the selected raw
% dataset plus participant-level XML files in the original open root.
metadata = struct();
paths = {};
datasetPath = scanPathValue(scan,'filelocation');
openRoot = '';
try
    if isobject(scan) && isprop(scan,'source_metadata') && isstruct(scan.source_metadata) && ...
            isfield(scan.source_metadata,'OpenRootPath')
        openRoot = textValue(scan.source_metadata.OpenRootPath);
    elseif isstruct(scan) && isfield(scan,'source_metadata') && isstruct(scan.source_metadata) && ...
            isfield(scan.source_metadata,'OpenRootPath')
        openRoot = textValue(scan.source_metadata.OpenRootPath);
    end
catch
end
if ~isempty(datasetPath)
    paths = [paths, {fullfile(datasetPath,'Serie.xml'), fullfile(datasetPath,'header.xml')}]; %#ok<AGROW>
end
if ~isempty(openRoot) && isfolder(openRoot)
    rootXml = dir(fullfile(openRoot,'*.xml'));
    for k = 1:numel(rootXml)
        paths{end+1} = fullfile(rootXml(k).folder,rootXml(k).name); %#ok<AGROW>
    end
    % Common participant/study metadata file names, including one level above
    % a sequence directory when direct-file mode was used.
    names = {'Patient.xml','patient.xml','Subject.xml','subject.xml', ...
        'Study.xml','study.xml','Session.xml','session.xml','Exam.xml','exam.xml'};
    for k = 1:numel(names)
        paths{end+1} = fullfile(openRoot,names{k}); %#ok<AGROW>
    end
end
if ~isempty(datasetPath)
    parentPath = fileparts(datasetPath);
    names = {'Patient.xml','patient.xml','Subject.xml','subject.xml', ...
        'Study.xml','study.xml','Session.xml','session.xml','Exam.xml','exam.xml'};
    for k = 1:numel(names)
        paths{end+1} = fullfile(parentPath,names{k}); %#ok<AGROW>
    end
end
paths = unique(paths,'stable');
for k = 1:numel(paths)
    fn = paths{k};
    if ~isfile(fn), continue; end
    try
        txt = fileread(fn);
        metadata = collectXmlIdentity(metadata,txt);
    catch
    end
end
end

function pathText = scanPathValue(scan,name)
pathText = '';
try
    if (isobject(scan) && isprop(scan,name)) || (isstruct(scan) && isfield(scan,name))
        pathText = textValue(scan.(name));
    end
catch
end
end

function metadata = collectXmlIdentity(metadata,txt)
% Parse key/value entries whether <value> has attributes or contains a
% parameter wrapper.  Header.xml and Serie.xml both use this pattern.
blocks = regexp(txt,'<entry(?:\s[^>]*)?>([\s\S]*?)</entry>','tokens');
for b = 1:numel(blocks)
    block = blocks{b}{1};
    kt = regexp(block,'<key(?:\s[^>]*)?>\s*([^<]+?)\s*</key>','tokens','once','ignorecase');
    if isempty(kt), continue; end
    vals = regexp(block,'<value(?:\s[^>]*)?>\s*([^<]+?)\s*</value>','tokens','ignorecase');
    value = '';
    for v = 1:numel(vals)
        candidate = xmlDecodeSimple(vals{v}{1});
        if ~isempty(candidate)
            value = candidate;
            break
        end
    end
    if ~isempty(value)
        metadata = assignIdentity(metadata,kt{1},value,false);
    end
end

% Direct patient/subject blocks sometimes store an internal numeric ID as a
% simple <id> element rather than as a metadata entry.
patientBlocks = regexp(txt,'<(patient|subject)(?:\s[^>]*)?>([\s\S]*?)</\1>','tokens','ignorecase');
for b = 1:numel(patientBlocks)
    block = patientBlocks{b}{2};
    leaf = regexp(block,'<([A-Za-z][A-Za-z0-9_-]*)(?:\s[^>]*)?>\s*([^<]+?)\s*</\1>', ...
        'tokens','ignorecase');
    firstName=''; lastName='';
    for q = 1:numel(leaf)
        key = leaf{q}{1}; value = xmlDecodeSimple(leaf{q}{2});
        nk = normalizeFieldName(key);
        if strcmp(nk,'FIRSTNAME'), firstName=value; end
        if strcmp(nk,'LASTNAME'), lastName=value; end
        metadata = assignIdentity(metadata,key,value,true);
    end
    if ~isfield(metadata,'PatientName') || isempty(metadata.PatientName)
        fullName = strtrim([lastName ' ' firstName]);
        if ~isempty(fullName), metadata.PatientName = fullName; end
    end
end

% Also catch explicit demographic tags outside a patient wrapper.
tagNames = {'PatientID','PatientIdentifier','PatientNumber','SubjectID','SubjectNumber', ...
    'PatientBirthDate','PatientBirthday','BirthDate','DateOfBirth','DOB', ...
    'PatientSex','SubjectSex','Sex','Gender','PatientWeight','SubjectWeight','WeightKg', ...
    'StudyID','StudyIdentifier','StudyNumber'};
for n = 1:numel(tagNames)
    pat = ['<' tagNames{n} '(?:\s[^>]*)?>\s*([^<]+?)\s*</' tagNames{n} '>'];
    tok = regexp(txt,pat,'tokens','once','ignorecase');
    if ~isempty(tok), metadata = assignIdentity(metadata,tagNames{n},xmlDecodeSimple(tok{1}),false); end
end
end

function metadata = assignIdentity(metadata,key,value,inPatientBlock)
if isempty(value), return; end
nk = normalizeFieldName(key);
field = '';
if any(strcmp(nk,{'PATIENT','PATIENTNAME','SUBJECTNAME','FULLNAME'}))
    field='PatientName';
elseif any(strcmp(nk,{'PATIENTID','PATIENTIDENTIFIER','PATIENTNUMBER','PATIENTINDEX', ...
        'SUBJECTID','SUBJECTIDENTIFIER','SUBJECTNUMBER','SUBJECTINDEX','PERSONID'})) || ...
        (inPatientBlock && strcmp(nk,'ID'))
    field='PatientID';
elseif strcmp(nk,'DOB') || strcmp(nk,'BIRTHDAY') || contains(nk,'BIRTHDATE') || contains(nk,'DATEOFBIRTH')
    field='PatientBirthDate';
elseif any(strcmp(nk,{'PATIENTSEX','SUBJECTSEX','SEX','PATIENTGENDER','SUBJECTGENDER','GENDER'}))
    field='PatientSex';
elseif any(strcmp(nk,{'PATIENTWEIGHT','SUBJECTWEIGHT','WEIGHTKG'})) || (inPatientBlock && strcmp(nk,'WEIGHT'))
    field='PatientWeight';
elseif any(strcmp(nk,{'STUDY','STUDYDESCRIPTION','STUDYNAME'}))
    field='StudyDescription';
elseif any(strcmp(nk,{'STUDYID','STUDYIDENTIFIER','STUDYNUMBER'}))
    field='StudyID';
elseif any(strcmp(nk,{'SESSION','SESSIONDATE'}))
    field='Session';
elseif any(strcmp(nk,{'PATIENTPOSITION','SUBJECTPOSITION'}))
    field='PatientPosition';
elseif any(strcmp(nk,{'SERIENAME','SERIESNAME','PROTOCOLNAME'}))
    field='SeriesDescription';
elseif any(strcmp(nk,{'SERIEID','SERIESID'}))
    field='SeriesID';
elseif any(strcmp(nk,{'ACQUISITIONDATE','ACQUISITIONDATETIME'}))
    field='AcquisitionDateTime';
end
if isempty(field), return; end
if ~isfield(metadata,field) || isempty(metadata.(field))
    metadata.(field) = value;
end
end

function out = xmlDecodeSimple(in)
out = strtrim(char(in));
out = regexprep(out,'<[^>]+>','');
out = strrep(out,'&amp;','&');
out = strrep(out,'&lt;','<');
out = strrep(out,'&gt;','>');
out = strrep(out,'&quot;','"');
out = strrep(out,'&apos;',char(39));
out = strtrim(out);
end

function out = normalizeKey(value)
out = lower(regexprep(strtrim(char(string(value))),'[^A-Za-z0-9]+',''));
end

function values = nonEmpty(values)
values = values(~cellfun(@isempty,values));
end
