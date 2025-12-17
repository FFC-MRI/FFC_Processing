function [obj] = correct_orientation(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

off1 = obj.param.OFF_CENTER_FIELD_OF_VIEW_1D;  % along 1D axis
off2 = obj.param.OFF_CENTER_FIELD_OF_VIEW_2D;  % along 2D axis
off3 = obj.param.OFF_CENTER_FIELD_OF_VIEW_3D;  % along 3D axis

FOV1 = obj.param.FIELD_OF_VIEW;
FOV2 = obj.param.FIELD_OF_VIEW_PHASE;
FOV3 = obj.param.FIELD_OF_VIEW_3D;


Nx = size(obj.complexkspace,1);
Ny = size(obj.complexkspace,2);
Nz = size(obj.complexkspace,3);

kx = ((0:Nx-1) - Nx/2) / FOV1;   % cycles/m along 1D axis
ky = ((0:Ny-1) - Ny/2) / FOV2;   % cycles/m along 2D axis
kz = ((0:Nz-1) - Nz/2) / FOV3;   % cycles/m along 3D axis

[KX, KY, KZ] = ndgrid(kx, ky, kz);  % matches K(i,j,k)

phase = exp(-1i * 2*pi * (KX*off1 + KY*off2 + KZ*0));
obj.complexkspace = obj.complexkspace .* phase;


     obj.complexkspace = permute(obj.complexkspace,[2 1 3 4 5 6 7]);

  obj.complexkspace = flipud(obj.complexkspace);

end

