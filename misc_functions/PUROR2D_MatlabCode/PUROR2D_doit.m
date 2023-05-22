%
function [out_x out_y] = PUROR2D_doit(phase_original,mask_2D,mask_2D_h,debug_PUROR)
%----------------------------------------------------------------------
    phase_range = 4.0*pi; %[-phase_range phase_range]
%----------------------------------------------------
%
[cell_signal_x cell_connect_x cell_seg_x] = dePULM_2D_ini(mask_2D);
[cell_signal_x_h cell_connect_x_h cell_seg_x_h] = dePULM_2D_ini(mask_2D_h);
%------------------------------------
[cell_signal_y cell_connect_y cell_seg_y] = dePULM_2D_ini_y(mask_2D); 
[cell_signal_y_h cell_connect_y_h cell_seg_y_h] = dePULM_2D_ini_y(mask_2D_h);
%-------------------------------------------------------------------
[phase_itoh_x phase_itoh_y] = dePULM_1D_itoh(phase_original);
%figure(10)
%subplot(4,2,1);imagesc(phase_itoh_x,[-phase_range phase_range]);colormap gray;axis square;axis off;
%subplot(4,2,2);imagesc(phase_itoh_y,[-phase_range phase_range]);colormap gray;axis square;axis off;
if debug_PUROR == 1
    figure(1)
    imshow(phase_itoh_x,[-phase_range phase_range]);
    figure(2)
    imshow(phase_itoh_y,[-phase_range phase_range]);
    figure(9)
    imshow((phase_itoh_y - phase_itoh_x),[-phase_range phase_range]);
    pause
end
%-------------------------------------------------------------------
[ori_phase_x] = dePULM_2D_mean(phase_itoh_x,cell_signal_x, cell_connect_x_h,cell_seg_x);
%----------------------------------------------
[ori_phase_y] = dePULM_2D_mean_y(phase_itoh_y,cell_signal_y, cell_connect_y_h,cell_seg_y);
%--------------------------------------------------------------------------
if debug_PUROR == 1   
    figure(3)
    imshow(ori_phase_x,[-phase_range phase_range]);
    figure(4)
    imshow(ori_phase_y,[-phase_range phase_range]);
    figure(9)
    imshow((ori_phase_y - ori_phase_x),[-phase_range phase_range]);
    pause
end    
%--------------------------------------------------------------------------    
%subplot(4,2,3);imagesc(ori_phase_x,[-phase_range phase_range]);colormap gray;axis square;axis off;
%subplot(4,2,4);imagesc(ori_phase_y,[-phase_range phase_range]);colormap gray;axis square;axis off;
%
    [q_br_y xy_start_dw xy_start_up] = dePULM_2D_quality(ori_phase_x, mask_2D_h);
%-----------------------------
    [unwrapped_phase_y] = dePULM_2D_xy_cen(ori_phase_x,ori_phase_y,xy_start_dw,xy_start_up,mask_2D_h);      
%-----------------------------------------------------------------
    [q_br_x xy_start_L xy_start_R] = dePULM_2D_quality_y(unwrapped_phase_y, mask_2D_h);
     if (xy_start_R - xy_start_L)/(xy_start_up - xy_start_dw) >= 0.5
    [unwrapped_phase_x] = dePULM_2D_xy_cen_y(unwrapped_phase_y,ori_phase_x,xy_start_L,xy_start_R,mask_2D_h);
     else
    [unwrapped_phase_x] = dePULM_refY_trlX_tmp(unwrapped_phase_y, ori_phase_x,mask_2D_h,xy_start_dw,xy_start_up);
     end
%-------------
%subplot(4,2,5);imagesc(unwrapped_phase_x,[-phase_range phase_range]);colormap gray;axis square;axis off;
%subplot(4,2,6);imagesc(unwrapped_phase_y,[-phase_range phase_range]);colormap gray;axis square;axis off;
if debug_PUROR == 1
    figure(5)
    imshow(unwrapped_phase_x,[-phase_range phase_range]);
    figure(6)
    imshow(unwrapped_phase_y,[-phase_range phase_range]);
    figure(9)
    imshow((unwrapped_phase_y - unwrapped_phase_x),[-phase_range phase_range]);
end    
%--------------------------------------------------------------------------
    seg_phi = pi;
    for test_loop = 1:6
    [unwrapped_phase_y] = dePULM_refX_trlY_final(unwrapped_phase_x, unwrapped_phase_y,seg_phi,cell_seg_y,cell_signal_y);
    %
    [unwrapped_phase_x] = dePULM_refY_trlX_final(unwrapped_phase_y, unwrapped_phase_x,seg_phi,cell_seg_x,cell_signal_x);
    %
    seg_phi = seg_phi/2;
        if debug_PUROR == 1 
        phase_range = 2.0*pi;
        figure(3)
        imshow(unwrapped_phase_x,[-phase_range phase_range]);
        figure(4)
        imshow(unwrapped_phase_y,[-phase_range phase_range]);
        figure(9)
        phase_range = 1.0*pi;
        imshow((unwrapped_phase_y - unwrapped_phase_x),[-phase_range phase_range]);
        test_loop
        pause
        end
    end
%-------------------------------------------------------------------------
    [out_y] = dePULM_2D_diff_y(unwrapped_phase_y);
    [out_x] = dePULM_2D_diff(unwrapped_phase_x); 
%   subplot(4,2,7);imagesc(out_x,[-phase_range phase_range]);colormap gray;axis square;axis off;
%   subplot(4,2,8);imagesc(out_y,[-phase_range phase_range]);colormap gray;axis square;axis off;
%--------------------------------------------------------------------   
%   phase_unwrapping was done
%--------------------------------------------------------------------------

