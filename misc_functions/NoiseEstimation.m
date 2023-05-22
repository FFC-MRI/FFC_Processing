%
function [th4noise] = NoiseEstimation(complex_image)
%
%-------------------------
   matrix_size = size(complex_image);
   if length(matrix_size) == 2
   matrix_size(3) = 1;
   end
    %
   for index_slice = 1:matrix_size(3)
   mag_original(:,:) = abs(complex_image(:,:,index_slice));  
   mag_vec = mag_original(:);
   mag_vec(isnan(mag_vec)) = [];
   mag_vec_gray = mat2gray(mag_vec);
   mag_vec(mag_vec_gray > 0.05) = [];
   mag_vec(mag_vec == 0) = [];
   th4noise(index_slice) = std(mag_vec);
   end
   %----------------------------------------------------------------------
