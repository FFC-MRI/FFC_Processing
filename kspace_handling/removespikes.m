function [kspace] = removespikes(kspace)
%Try and remove dubious kspace lines
%   Detailed explanation goes here

% if kspace(:,1,1) == 0
%     kspace(:,1,1) = kspace(:,2,1); 
% end
dims = size(kspace);
kspace = reshape(kspace,dims(1),dims(2),[]);
% for p=2:size(kspace,3)
%     for n=1:size(kspace,2)
%         thresh = std(kspace(:,:,:));
%         thresh = mean(squeeze(thresh(:,n,1:2:end)));
%         if std(kspace(:,n,p))>thresh*5
%             if p ==1
%                 kspace(:,n,p) =  kspace(:,n,p+1);
%             else
%                 kspace(:,n,p) =  kspace(:,n,p-1);
%             end
%         end
%     end
% end
if size(kspace,3) > 1

    kspacer = filloutliers(real(kspace),'spline','median',3,'ThresholdFactor',100);
    kspacei = filloutliers(imag(kspace),'spline','median',3,'ThresholdFactor',100);
    kspace = kspacer +1i*kspacei;


end

kspace = reshape(kspace,dims);

end

