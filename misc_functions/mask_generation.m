%
function [mask4unwrap, mask4supp, mask4stack] = mask_generation(complex_image,th4noise,th4unwrap,th4supp,th4stack)
%-------------------------
matrix_size = size(complex_image);
if length(matrix_size) == 2
    matrix_size(3) = 1;
end
%-------------------------
   H = fspecial('gaussian',3);  
%-------------------------   
   mask4unwrap = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
   mask4supp = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
   mask4stack = zeros(matrix_size(1),matrix_size(2),matrix_size(3));
   for index_slice = 1:matrix_size(3)
   complex_tmp(:,:) = complex_image(:,:,index_slice);  
   mag_original = abs(complex_tmp);
   background_noise = th4noise(index_slice);
   %--------------------------------------
   mask_tmp = ones(matrix_size(1),matrix_size(2));
   mask_tmp(mag_original < background_noise*th4stack) = 0;
   [L_0 num_region] = bwlabel(mask_tmp,4);
   mask_tmp(L_0 == 0) = 0;
   %figure(3)
   %imagesc(mask_tmp);colormap gray;axis square;axis off;
   mask4stack(:,:,index_slice) = mask_tmp;     
   %--------------------------------------
   complex_tmp = imfilter(complex_tmp,H);
   complex_tmp(isnan(complex_tmp)) = 0;
   mag_original = abs(complex_tmp);
   mask_tmp_0 = ones(matrix_size(1),matrix_size(2));
   mask_tmp_0(mag_original < background_noise*th4unwrap) = 0;
   [L_0] = bwlabel(mask_tmp_0,4);
   mask_tmp_0(L_0 == 0) = 0;
   %
   mask4unwrap(:,:,index_slice) = mask_tmp_0;
   %--------------------------------------------
   mask_tmp = ones(matrix_size(1),matrix_size(2));
   mask_tmp(mag_original < background_noise*th4supp) = 0;       
   mask_tmp(mask_tmp_0 == 0) = 0;
   L_h = bwlabel(mask_tmp,4);
   mask_tmp(L_h == 0) = 0;
   mask4supp(:,:,index_slice) = mask_tmp;    
   %
   %index_slice
   %pause
   end
   %---------------------------------------------------------------------
