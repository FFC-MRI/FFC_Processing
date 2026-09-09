%% Denoise the image
image = squeeze(image);
imagen = image./max(image(:)); %normalise images

N_times = size(image,3);
N_fields = size(image,4);
figure;imagesc(imagen(:,:,1,1)); %plot image to select the noise area

     for i = 1:N_fields
         for  j = 1:N_times
             img = imagen(:,:,j,i);
             noise = std2(img(1:20,1:100))%+ 0.2*std2(img(1:15,1:100)); %(y,x)
             denoised_image = BM3D(img,noise);
             denoised_image_5D(:,:,1,j,i) = denoised_image;
        end
     end
     
denoised_image = denoised_image_5D;

%%  %% plot images
figure
imagesc(abs(image(:,:,1,1)));colormap(gray);%camroll(90);
title('Raw')

figure
imagesc(abs(img(:,:,1,1,1)));colormap(gray);%camroll(90);
title('Denoised')


image = squeeze(image);
figure;
i = 1;
for nBevo = 1: size(image,4)
    for nTevo = 1:size(image,3)
        subplot(size(image,4),size(image,3),i)
        imagesc(abs(image(:,:,nTevo,nBevo)));
        i = i+1;
        axis off
        
    end
  end
colormap(gray);


denoised_image = squeeze(denoised_image);
figure;

i = 1;
for nBevo = 1: size(denoised_image,4)
    for nTevo = 1:size(denoised_image,3)
        subplot(size(denoised_image,4),size(denoised_image,3),i)
        imagesc(abs(denoised_image(:,:,nTevo,nBevo)));
        i = i+1;
        axis off
    
    end
end
colormap(gray);


