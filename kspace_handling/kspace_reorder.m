function [A] = kspace_reorder(kspace,obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sequence = obj.sequence;




switch sequence
    case "h9_flash_v2"
        table = [0,-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8,8,-9,9,-10,10,-11,11,-12,12,-13,13,-14,14,-15,15,-16,16,-17,17,-18,18,-19,19,-20,20,-21,21,-22,22,-23,23,-24,24,-25,25,-26,26,-27,27,-28,28,-29,29,-30,30,-31,31,-32,32,-33,33,-34,34,-35,35,-36,36,-37,37,-38,38,-39,39,-40,40,-41,41,-42,42,-43,43,-44,44,-45,45,-46,46,-47,47,-48,48,-49,49,-50,50,-51,51,-52,52,-53,53,-54,54,-55,55,-56,56,-57,57,-58,58,-59,59,-60,60,-61,61,-62,62,-63,63,-64];
                       tempdata = kspace;
                       
                       tempmatrix = zeros(size(tempdata));
                       for n=1:size(tempdata,2)  
                       tempmatrix(:,1+size(tempdata,2)/2+table(n),:,:,:,:,:,:,:) =  tempdata(:,n,:,:,:,:,:,:,:);
                       end
                       
            A = tempmatrix;
            
            case "H9_ir_se"
                kspace(:,:,1,:,:,:,:,:) = kspace(:,:,2,:,:,:,:,:);
                A = kspace;
    otherwise
        
   
            
                
        A= kspace;
end

end

