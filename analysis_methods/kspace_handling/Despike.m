function [kspace] = Despike(kspace)
%Try and remove dubious kspace lines
%   Detailed explanation goes here

% if kspace(:,1,1) == 0
%     kspace(:,1,1) = kspace(:,2,1); 
% end
dims = size(kspace);
kspace = reshape(kspace,dims(1),dims(2),[]);
for p=1:size(kspace,3)
kspace = despike_kspace_sym(kspace(:,:,p),5);
end



end

kspace = reshape(kspace,dims);



