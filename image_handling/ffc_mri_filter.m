function [filtered_images] = ffc_mri_filter(images,filter_type,kernel)
%P
%   Detailed explanation goes here


% % modification:  normalising the images field-by-field (LB 1/03/20, needs
% % optimising)
% for b = 1:size(images,5)
%     images(:,:,:,:,b) = images(:,:,:,:,b)./repmat(images(:,:,:,1,1),1,1,1,size(images,4));
% end

dim = size(images);
tempimages = abs(reshape(images,dim(1),dim(2),[])); %reshape for processing ease

switch filter_type
    case 'Generalised Total Variation'
        parfor n=1:size(tempimages,3)
            filtered_images(:,:,n) = imtgvsmooth(tempimages(:,:,n),kernel,kernel,200);
        end
        
    case 'Total Variation'
        parfor n=1:size(tempimages,3)
            filtered_images(:,:,n) = TVL1denoise(tempimages(:,:,n),kernel,100);
        end
    case 'Deep Learning'
        net = denoisingNetwork('DnCNN');
        parfor n=1:size(tempimages,3)
            tempimages(:,:,n) =  tempimages(:,:,n)./max(max( tempimages(:,:,n)));
              filtered_images(:,:,n) = denoiseImage(tempimages(:,:,n), net);
        end
    otherwise
        filtered_images = tempimages;
end
filtered_images = reshape(filtered_images,dim);
end

