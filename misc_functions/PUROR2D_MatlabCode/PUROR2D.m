function [ unwrapped_phase ] = PUROR2D( complex_image,mask4unwrap, mask4supp, mask4stack, debug_PUROR )
%UNTITLED Summary of this function goes here
%Liu J and Drangova M, MRM, Vol 68(4): 1303-1316, 2012
matrix_size = size(complex_image);
if length(matrix_size) == 2
matrix_size(3) = 1;    
end
%
if debug_PUROR == 1
display_range = 4.0*pi;
end
%
%----------------------------------------------------------------------
%   Start phase unwrapping 2D
%----------------------------------------------------------------------
phase_original = zeros(matrix_size);
unwrapped_phase_x = zeros(matrix_size);
unwrapped_phase_y = zeros(matrix_size);
for index_slice = 1:matrix_size(3)    
    data_original(:,:) = complex_image(:,:,index_slice); 
    phase_slice = angle(data_original);
    phase_slice(isnan(phase_slice)) = 0;
    phase_original(:,:,index_slice) = phase_slice;
    %----------------
    [unwrapped_x unwrapped_y] = PUROR2D_doit(phase_slice,mask4unwrap(:,:,index_slice),mask4supp(:,:,index_slice),debug_PUROR);
    unwrapped_phase_x(:,:,index_slice) = unwrapped_x;
    unwrapped_phase_y(:,:,index_slice) = unwrapped_y;    
end
%--------------------------------------------------------------------------
if debug_PUROR == 1
    for index_slice = 1:matrix_size(3)
        figure(1)
        subplot(1,3,1);imagesc(phase_original(:,:,index_slice),[-display_range display_range]);colormap gray;axis square;axis off;
        subplot(1,3,2);imagesc(unwrapped_phase_x(:,:,index_slice),[-display_range display_range]);colormap gray;axis square;axis off;
        subplot(1,3,3);imagesc(unwrapped_phase_y(:,:,index_slice),[-display_range display_range]);colormap gray;axis square;axis off;
        index_slice
        pause
    end
end
%--------------------------------------------------------------------------
%   2D phase_unwrapping was done
%--------------------------------------------------------------------------
%stacking the multiple slices
%--------------------------------------------------------------------------
if matrix_size(3) > 1
[ unwrapped_phase_x, unwrapped_phase_y ] = StackSlice(unwrapped_phase_x, unwrapped_phase_y,mask4unwrap, mask4stack);    
end
%---------------------------------------------------------------------
if debug_PUROR == 1
    for index_slice = 1:matrix_size(3)
        figure(1)
        subplot(1,3,1);imagesc(phase_original(:,:,index_slice),[-display_range display_range]);colormap gray;axis square;axis off;
        subplot(1,3,2);imagesc(unwrapped_phase_x(:,:,index_slice),[-display_range display_range]);colormap gray;axis square;axis off;
        subplot(1,3,3);imagesc(unwrapped_phase_y(:,:,index_slice),[-display_range display_range]);colormap gray;axis square;axis off;
        index_slice
        pause
    end
end
%---------------
unwrapped_phase = unwrapped_phase_x;
%---------------
end

