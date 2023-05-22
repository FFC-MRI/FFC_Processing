function [ unwrapped_phase_x, unwrapped_phase_y ] = StackSlice(unwrapped_phase_x, unwrapped_phase_y,mask4unwrap, mask4stack)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pi_2 = 2*pi;
matrix_size = size(unwrapped_phase_x);
for index_slice = 1:matrix_size(3)
MASK_2D_ori(:,:) = mask4unwrap(:,:,index_slice);
MASK_2D_ori_h(:,:) = mask4stack(:,:,index_slice);
MASK_stas(:,:) = mask4stack(:,:,index_slice);
%
phase_tmp_mean = unwrapped_phase_x(:,:,index_slice);
phase_tmp_mean = phase_tmp_mean(:);
phase_tmp_mean(isnan(phase_tmp_mean)) = 0;
phase_tmp_mean(MASK_stas(:) == 0) = [];
phase_tmp_mean(isnan(phase_tmp_mean)) = [];
tmp_mean = mean(phase_tmp_mean);
%tmp_mean/pi
unwrapped_phase_x = unwrapped_phase_x - pi_2*round(tmp_mean/pi_2);
%
tmp_mean = tmp_mean - pi_2*round(tmp_mean/pi_2);
mean_2D_x(index_slice) = tmp_mean;
%
phase_tmp_mean =  unwrapped_phase_y(:,:,index_slice);
phase_tmp_mean = phase_tmp_mean(:);
phase_tmp_mean(isnan(phase_tmp_mean)) = 0;
phase_tmp_mean(MASK_stas(:) == 0) = [];
phase_tmp_mean(isnan(phase_tmp_mean)) = [];
tmp_mean = mean(phase_tmp_mean);
%tmp_mean/pi
unwrapped_phase_y = unwrapped_phase_y - pi_2*round(tmp_mean/pi_2);
tmp_mean = tmp_mean - pi_2*round(tmp_mean/pi_2);
mean_2D_y(index_slice) = tmp_mean;
    %   
end % for index_slice
%----------------------------------
slice_centre = round(matrix_size(3)/2);
unwrap_mean_x = unwrap(mean_2D_x);
unwrap_mean_x = unwrap_mean_x - pi_2*round(unwrap_mean_x(slice_centre)/pi_2);
diff_mean_x = mean_2D_x - unwrap_mean_x;
for index_slice = 1:matrix_size(3)
unwrapped_phase_x(:,:,index_slice) = unwrapped_phase_x(:,:,index_slice) - pi_2*round(diff_mean_x(index_slice)/pi_2);
end
%---------------
unwrap_mean_y = unwrap(mean_2D_y);
unwrap_mean_y = unwrap_mean_y - pi_2*round(unwrap_mean_y(slice_centre)/pi_2);
diff_mean_y = mean_2D_y - unwrap_mean_y;
for index_slice = 1:matrix_size(3)
unwrapped_phase_y(:,:,index_slice) = unwrapped_phase_y(:,:,index_slice) - pi_2*round(diff_mean_y(index_slice)/pi_2);
end
%-----------------------------------------
end%


