function K = ffc_read_validation_fixture(datFile, shape)
%FFC_READ_VALIDATION_FIXTURE Read an RS2D Cartesian export for batch QA.
% Pass the dimensions verified from header.xml. Uses the app's [imag,real]
% big-endian single convention. Use the app loader for general data.
if nargin < 2, shape = [128 128 8 7 10]; end
fid = fopen(datFile, 'rb', 'ieee-be');
assert(fid >= 0, 'Unable to open RS2D data.dat');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
a = fread(fid, [2 prod(shape)], '*single');
assert(numel(a)==2*prod(shape) && isempty(fread(fid,1,'*uint8')), ...
    'RS2D data length disagrees with the supplied shape.');
K = reshape(complex(a(2,:),a(1,:)),shape);
end
