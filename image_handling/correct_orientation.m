function [obj] = correct_orientation(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if obj.phase_direction == 0
     obj.complexkspace = permute(obj.complexkspace,[2 1 3 4 5 6 7]);
end
  obj.complexkspace = flipud(obj.complexkspace);

end

