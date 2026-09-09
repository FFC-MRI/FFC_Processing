function cvi = cvi_read_workspace(filename)
%CVI_READ_WORKSPACE Read Circle cvi42 contour coordinates from a workspace.
%
%   CVI = CVI_READ_WORKSPACE(FILE) supports the zlib-compressed .cvi42ws
%   workspace format used by cvi42 5.x/early 6.x and the XML .cvi42wsx
%   export format.  The binary reader deliberately extracts only the stable
%   fields needed by FFC Processing: contour name, image key, sub-pixel
%   resolution and contour points.  It does not attempt to deserialize the
%   complete proprietary workspace object graph.
%
%   Returned contour point coordinates are in the source DICOM image pixel
%   coordinate system.  cvi42 uses zero-based sub-pixel coordinates, so no
%   MATLAB one-based offset is applied here.

    arguments
        filename (1,:) char
    end
    if ~isfile(filename)
        error('cvi_read_workspace:FileNotFound','Workspace not found: %s',filename);
    end

    [~,~,ext] = fileparts(filename);
    switch lower(ext)
        case '.cvi42ws'
            contours = readBinaryWorkspace(filename);
        case '.cvi42wsx'
            contours = readXmlWorkspace(filename);
        otherwise
            error('cvi_read_workspace:UnsupportedFormat', ...
                'Expected a .cvi42ws or .cvi42wsx file.');
    end

    if isempty(contours)
        error('cvi_read_workspace:NoContours', ...
            'No supported SAX ventricular contours were found in the workspace.');
    end

    cvi = struct();
    cvi.workspaceFile = filename;
    cvi.contours = contours;
    cvi.importedUtc = datetime('now','TimeZone','UTC');
end

function contours = readBinaryWorkspace(filename)
    fid = fopen(filename,'rb');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    packed = fread(fid,inf,'*uint8').';
    if isempty(packed)
        error('cvi_read_workspace:EmptyWorkspace','The workspace is empty.');
    end
    try
        bytes = cviZlibDecode(packed);
    catch err
        error('cvi_read_workspace:DecompressionFailed', ...
            ['Could not decompress the cvi42 workspace. The importer tried ', ...
             'MATLAB zlibdecode (when available) and a Java zlib fallback: %s'], ...
             err.message);
    end
    bytes = uint8(bytes(:).');

    supported = { ...
        'saendocardialContour','LVEndo'; ...
        'saepicardialContour','LVEpi'; ...
        'sarvendocardialContour','RVEndo'};
    uidPattern = utf16beAscii('1.3.');
    uidStarts = strfind(bytes,uidPattern);
    if isempty(uidStarts)
        error('cvi_read_workspace:NoImageKeys', ...
            'The workspace did not contain recognisable cvi42 image keys.');
    end

    template = struct('Type','','CVIName','','ImageKey','','ImageUID','', ...
        'FrameIndex0',NaN,'SubpixelResolution',NaN,'PointsPixel',[], ...
        'PixelArea',NaN,'WorkspaceOffset',NaN,'DicomFile','', ...
        'PhaseIndex',NaN,'PointsPatientMm',[],'SliceIndex',NaN, ...
        'DicomPhaseIndex',NaN,'FramePositionPatientMm',[NaN NaN NaN], ...
        'FrameNormalPatient',[NaN NaN NaN]);
    contours = repmat(template,0,1);

    for typeIndex = 1:size(supported,1)
        cviName = supported{typeIndex,1};
        typeName = supported{typeIndex,2};
        namePattern = utf16beAscii(cviName);
        starts = strfind(bytes,namePattern);
        for occurrence = 1:numel(starts)
            contourStart = starts(occurrence);
            preceding = uidStarts(uidStarts < contourStart);
            if isempty(preceding)
                continue;
            end
            imageKey = extractUtf16Ascii(bytes,preceding(end));
            if ~startsWith(imageKey,'1.3.')
                continue;
            end

            token = regexp(imageKey,'^(.*)multi(\d+)$','tokens','once');
            if isempty(token)
                imageUID = imageKey;
                frameIndex0 = 0;
            else
                imageUID = token{1};
                frameIndex0 = str2double(token{2});
            end

            afterName = contourStart + numel(namePattern);
            if afterName + 3 > numel(bytes)
                continue;
            end
            subpixel = double(readBEUInt32(bytes,afterName));
            if ~isfinite(subpixel) || subpixel < 1 || subpixel > 1024
                subpixel = 1;
            end

            searchLast = min(numel(bytes),afterName + 512);
            freeRel = strfind(bytes(afterName:searchLast),uint8([double('FREE') 0]));
            if isempty(freeRel)
                continue;
            end
            freeStart = afterName + freeRel(1) - 1;
            if freeStart + 12 > numel(bytes)
                continue;
            end
            payloadLength = double(readBEUInt32(bytes,freeStart + 5));
            pointCount = double(readBEUInt32(bytes,freeStart + 9));
            if pointCount < 3 || pointCount > 1e6 || ...
                    payloadLength ~= 4 + 8 .* pointCount
                continue;
            end
            dataStart = freeStart + 13;
            dataEnd = dataStart + 8 .* pointCount - 1;
            if dataEnd > numel(bytes)
                continue;
            end
            raw = typecast(uint8(bytes(dataStart:dataEnd)),'int32');
            raw = swapbytes(raw);
            raw = reshape(double(raw),2,[]).';
            points = raw ./ subpixel;

            item = template;
            item.Type = typeName;
            item.CVIName = cviName;
            item.ImageKey = imageKey;
            item.ImageUID = imageUID;
            item.FrameIndex0 = frameIndex0;
            item.SubpixelResolution = subpixel;
            item.PointsPixel = points;
            item.PixelArea = polyarea(points(:,1),points(:,2));
            item.WorkspaceOffset = contourStart;
            contours(end+1,1) = item; %#ok<AGROW>
        end
    end

    if ~isempty(contours)
        [~,order] = sort([contours.WorkspaceOffset]);
        contours = contours(order);
    end
end


function bytes = cviZlibDecode(packed)
%CVI_ZLIB_DECODE Decode a standard RFC1950/zlib stream on old and new MATLAB.
%
% Newer MATLAB releases provide zlibdecode directly. Older releases do not.
% The Java fallback deliberately avoids InputStream.read(MATLAB_ARRAY,...): on
% some MATLAB releases Java writes into a temporary copy of a MATLAB numeric
% array, leaving the MATLAB buffer unchanged (typically all zeros). That was
% the cause of TEST4 successfully "decompressing" but then finding no cvi42
% keys. All bulk-copy paths below now keep the copy inside Java; the final
% fallback uses InputStream.read() one byte at a time and is slower but is
% compatible with old MATLAB/JVM combinations.

    packed = uint8(packed(:));

    % Prefer MATLAB's implementation when available.
    if exist('zlibdecode','file') == 2 || exist('zlibdecode','builtin') == 5
        bytes = zlibdecode(packed);
        bytes = uint8(bytes(:));
        return;
    end

    if ~usejava('jvm')
        error('cvi_read_workspace:NoZlibDecoder', ...
            ['This MATLAB release does not provide zlibdecode and the JVM ', ...
             'is unavailable, so the .cvi42ws file cannot be decompressed.']);
    end

    % Java byte is signed. typecast preserves the original 8-bit pattern;
    % numerical uint8->int8 conversion would clip values above 127.
    signedInput = typecast(packed,'int8');

    % ---- Fast path 1: Java 9+ InputStream.readAllBytes -----------------
    % Everything is copied inside Java, so there is no MATLAB array pass-by-
    % value problem. Older MATLAB releases often ship Java 8, so this is
    % optional.
    try
        input = javaObject('java.io.ByteArrayInputStream', signedInput);
        inflater = javaObject('java.util.zip.InflaterInputStream', input);
        cleanup = onCleanup(@() closeJavaInputStreams(inflater,input)); %#ok<NASGU>
        signedOutput = inflater.readAllBytes();
        bytes = javaByteArrayToUint8(signedOutput);
        if ~isempty(bytes)
            return;
        end
    catch
        % Re-open a fresh stream for the next method.
    end

    % ---- Fast path 2: Apache Commons IO, bundled with many MATLABs ------
    try
        input = javaObject('java.io.ByteArrayInputStream', signedInput);
        inflater = javaObject('java.util.zip.InflaterInputStream', input);
        cleanup = onCleanup(@() closeJavaInputStreams(inflater,input)); %#ok<NASGU>
        signedOutput = javaMethod('toByteArray','org.apache.commons.io.IOUtils',inflater);
        bytes = javaByteArrayToUint8(signedOutput);
        if ~isempty(bytes)
            return;
        end
    catch
        % Commons IO is not guaranteed to be on every MATLAB Java classpath.
    end

    % ---- Compatibility path: scalar Java reads -------------------------
    % This is intentionally simple and reliable. A typical cvi42 workspace
    % expands to around 1 MB, so the one-off import cost is acceptable.
    input = javaObject('java.io.ByteArrayInputStream', signedInput);
    inflater = javaObject('java.util.zip.InflaterInputStream', input);
    cleanup = onCleanup(@() closeJavaInputStreams(inflater,input)); %#ok<NASGU>

    capacity = max(1048576, 4*numel(packed));
    bytes = zeros(capacity,1,'uint8');
    n = 0;
    while true
        value = inflater.read();
        if value < 0
            break;
        end
        n = n + 1;
        if n > numel(bytes)
            bytes(end+1:end+max(1048576,numel(bytes))) = uint8(0); %#ok<AGROW>
        end
        bytes(n) = uint8(value);
    end
    bytes = bytes(1:n);

    if isempty(bytes)
        error('cvi_read_workspace:JavaZlibEmpty', ...
            'The Java zlib decoder returned an empty stream.');
    end
end

function bytes = javaByteArrayToUint8(javaBytes)
%JAVA_BYTE_ARRAY_TO_UINT8 Preserve the bit pattern of Java signed bytes.
    signed = int8(javaBytes);
    bytes = typecast(signed(:),'uint8');
end

function closeJavaInputStreams(inflater,input)
    try, inflater.close(); catch, end
    try, input.close(); catch, end
end

function contours = readXmlWorkspace(filename)
    % cvi42wsx is an XML hash.  This implementation follows the public
    % cvi42wsx convention used by research tools such as shortCardiac.
    doc = xmlread(filename);
    template = struct('Type','','CVIName','','ImageKey','','ImageUID','', ...
        'FrameIndex0',NaN,'SubpixelResolution',NaN,'PointsPixel',[], ...
        'PixelArea',NaN,'WorkspaceOffset',NaN,'DicomFile','', ...
        'PhaseIndex',NaN,'PointsPatientMm',[],'SliceIndex',NaN, ...
        'DicomPhaseIndex',NaN,'FramePositionPatientMm',[NaN NaN NaN], ...
        'FrameNormalPatient',[NaN NaN NaN]);
    contours = repmat(template,0,1);
    imageStates = findHashKeyNodes(doc.getDocumentElement(),'ImageStates');
    supported = containers.Map( ...
        {'saendocardialContour','saepicardialContour','sarvendocardialContour'}, ...
        {'LVEndo','LVEpi','RVEndo'});
    for n = 1:numel(imageStates)
        imageStateContainer = imageStates{n};
        imageNodes = elementChildren(imageStateContainer);
        for i = 1:numel(imageNodes)
            imageKey = char(imageNodes{i}.getAttribute('Hash:key'));
            contourContainers = findHashKeyNodes(imageNodes{i},'Contours');
            for cc = 1:numel(contourContainers)
                contourNodes = elementChildren(contourContainers{cc});
                for k = 1:numel(contourNodes)
                    cviName = char(contourNodes{k}.getAttribute('Hash:key'));
                    if ~isKey(supported,cviName), continue; end
                    pointNodes = findHashKeyNodes(contourNodes{k},'Points');
                    subNodes = findHashKeyNodes(contourNodes{k},'SubpixelResolution');
                    if isempty(pointNodes), continue; end
                    subpixel = 1;
                    if ~isempty(subNodes)
                        subpixel = str2double(char(subNodes{1}.getTextContent()));
                        if ~isfinite(subpixel) || subpixel <= 0, subpixel = 1; end
                    end
                    children = elementChildren(pointNodes{1});
                    points = zeros(numel(children),2);
                    keep = false(numel(children),1);
                    for p = 1:numel(children)
                        xs = children{p}.getElementsByTagName('Point:x');
                        ys = children{p}.getElementsByTagName('Point:y');
                        if xs.getLength == 0 || ys.getLength == 0, continue; end
                        points(p,1) = str2double(char(xs.item(0).getTextContent()));
                        points(p,2) = str2double(char(ys.item(0).getTextContent()));
                        keep(p) = all(isfinite(points(p,:)));
                    end
                    points = points(keep,:) ./ subpixel;
                    if size(points,1) < 3, continue; end
                    token = regexp(imageKey,'^(.*)multi(\d+)$','tokens','once');
                    if isempty(token)
                        imageUID = imageKey; frameIndex0 = 0;
                    else
                        imageUID = token{1}; frameIndex0 = str2double(token{2});
                    end
                    item = template;
                    item.Type = supported(cviName);
                    item.CVIName = cviName;
                    item.ImageKey = imageKey;
                    item.ImageUID = imageUID;
                    item.FrameIndex0 = frameIndex0;
                    item.SubpixelResolution = subpixel;
                    item.PointsPixel = points;
                    item.PixelArea = polyarea(points(:,1),points(:,2));
                    contours(end+1,1) = item; %#ok<AGROW>
                end
            end
        end
    end
end

function nodes = findHashKeyNodes(node,key)
    nodes = {};
    if isempty(node), return; end
    if node.getNodeType == node.ELEMENT_NODE && ...
            strcmp(char(node.getAttribute('Hash:key')),key)
        nodes{end+1} = node; %#ok<AGROW>
    end
    children = node.getChildNodes();
    for k = 0:children.getLength()-1
        child = children.item(k);
        if child.getNodeType == child.ELEMENT_NODE
            found = findHashKeyNodes(child,key);
            nodes = [nodes,found]; %#ok<AGROW>
        end
    end
end

function children = elementChildren(node)
    children = {};
    nodes = node.getChildNodes();
    for k = 0:nodes.getLength()-1
        child = nodes.item(k);
        if child.getNodeType == child.ELEMENT_NODE
            children{end+1} = child; %#ok<AGROW>
        end
    end
end

function pattern = utf16beAscii(text)
    ascii = uint8(char(text));
    pattern = reshape([zeros(1,numel(ascii),'uint8'); ascii],1,[]);
end

function text = extractUtf16Ascii(bytes,startIndex)
    chars = char.empty(1,0);
    index = startIndex;
    while index + 1 <= numel(bytes) && bytes(index) == 0 && ...
            bytes(index+1) >= 32 && bytes(index+1) <= 126
        chars(end+1) = char(bytes(index+1)); %#ok<AGROW>
        index = index + 2;
    end
    text = chars;
end

function value = readBEUInt32(bytes,startIndex)
    b = uint32(bytes(startIndex:startIndex+3));
    value = bitshift(b(1),24) + bitshift(b(2),16) + ...
        bitshift(b(3),8) + b(4);
end
