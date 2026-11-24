function uid = ffc_make_dicom_uid(key)
%FFC_MAKE_DICOM_UID Create a DICOM UID, deterministically when KEY is supplied.
% Uses the UUID-derived 2.25 root.  Stable keys are useful for keeping a
% participant/study FrameOfReferenceUID consistent across separately exported
% FCI sequences without consulting any external DICOM files.

if nargin < 1, key = ''; end
key = char(string(key));
try
    if ~isempty(key)
        md = java.security.MessageDigest.getInstance('SHA-256');
        raw = md.digest(uint8(unicode2native(key,'UTF-8')));
        bytes = mod(double(raw(:)),256);
        hex = lower(reshape(dec2hex(bytes(1:16),2).',1,[]));
        big = javaObject('java.math.BigInteger',hex,16);
        uid = ['2.25.' char(big.toString())];
        return
    end
catch
end
try
    if isempty(key) && exist('dicomuid','file') == 2
        uid = dicomuid;
        return
    end
catch
end
try
    jUUID = javaMethod('randomUUID','java.util.UUID');
    hex = strrep(char(jUUID.toString()),'-','');
    big = javaObject('java.math.BigInteger',hex,16);
    uid = ['2.25.' char(big.toString())];
catch
    c = clock;
    seed = abs(round(sum(c.*[1e10 1e8 1e6 1e4 1e2 1])));
    uid = ['2.25.' sprintf('%.0f',seed) sprintf('%09d',randi(1e9-1))];
end
end
