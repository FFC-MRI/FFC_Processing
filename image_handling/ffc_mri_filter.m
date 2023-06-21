function [filtered_images] = ffc_mri_filter(images,filter_type,kernel)
% This function decides which noise filter to use based on the input
% 'filter_type' from the GUI drop down.


% % modification:  normalising the images field-by-field (LB 1/03/20, needs
% % optimising)
% for b = 1:size(images,5)
%     images(:,:,:,:,b) = images(:,:,:,:,b)./repmat(images(:,:,:,1,1),1,1,1,size(images,4));
% end

previousLocation = pwd;
[localDirectory,~,~] = fileparts(mfilename('fullpath'));
cd(localDirectory);
cd ..\image_handling\noise_filters;
Files = dir('*.m');
num_files = length(Files);
filterArray = {};

for i=1:num_files
    [pathstr, name, ext] = fileparts(Files(i).name);
    if strcmp(name,'NoiseFilterTemplate') == false
        filterArray{end+1} = feval(name);
    end
end
cd(previousLocation);


dim = size(images);
tempimages = abs(reshape(images,dim(1),dim(2),[])); %reshape for processing ease

filtered_images = tempimages;


for i=1:length(filterArray)

    if strcmp(filter_type,filterArray{i}.processName)== true
        filtered_images = filterArray{i}.Filter(tempimages,kernel);
    end
end


filtered_images = reshape(filtered_images,dim);
end

