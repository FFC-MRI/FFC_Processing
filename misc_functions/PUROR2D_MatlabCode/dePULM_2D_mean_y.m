%
function [out_2D] = dePULM_2D_mean_y(phase_original,cell_signal, cell_connect,cell_seg_y)
%----------------------------------------------------------------------
    pi_2 = pi*2.0;
    [yr_tmp xr_tmp] = size(phase_original);
    cen_line = round(xr_tmp/2);
%----------------------------------------------------------------
    % generate the mask-out unwrapped_phase_y
    phase_tmp = zeros(yr_tmp,xr_tmp);
%--------------------------------------------------------------------
    for index_x = 1:xr_tmp
        index_u = cell_signal{1,index_x};
        len_u = length(index_u);
        if len_u > yr_tmp/4
        phase_u = phase_original(index_u,index_x);       
        g_seg = cell_seg_y{1,index_x};
        phi_good = 0.5*pi;
        [index_ls] = dePULM_1D_ls(phase_u,phi_good,g_seg);
        [phase_tmp_1] = dePULM_1D(phase_u,index_u,index_ls);
        phase_tmp(index_u,index_x) = phase_tmp_1(1:len_u);
        else
        phase_tmp(index_u,index_x) = phase_original(index_u,index_x);            
        end      
        %-------
    end    
    %------------------------------------------------------------------
    mean_connect = zeros(1,xr_tmp);
    for index_x = 2:xr_tmp %
        phase_y = phase_tmp(:,index_x);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   calculating the connect_mean
        index_c = cell_connect{1,index_x};
        index_s = cell_signal{1,index_x};
        if isempty(index_s)
        mean_connect(index_x) = mean_connect(index_x - 1);        
        else
        %
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_x) = mean_liu(phase_y(index_s));        
            else
            mean_connect(index_x) = mean_liu(phase_y(index_c));
            end
        %
            phase_tmp(index_s,index_x) = phase_tmp(index_s,index_x) - round(mean_connect(index_x)/pi_2)*pi_2;
        %
            phase_y = phase_tmp(:,index_x);
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_x) = mean_liu(phase_y(index_s));        
            else
            mean_connect(index_x) = mean_liu(phase_y(index_c));
            end
        %
        end
    end
    %-------------------------------------------------------- 
    index_x = 1;
        phase_y = phase_tmp(:,index_x);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   calculating the connect_mean
        index_c = cell_connect{1,index_x};
        index_s = cell_signal{1,index_x};
        if isempty(index_s)
        mean_connect(index_x) = mean_connect(index_x + 1);        
        else
        %
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_x) = mean_liu(phase_y(index_s));        
            else
            mean_connect(index_x) = mean_liu(phase_y(index_c));
            end
        %
            phase_tmp(index_s,index_x) = phase_tmp(index_s,index_x) - round(mean_connect(index_x)/pi_2)*pi_2;
        %
            phase_y = phase_tmp(:,index_x);
            if isempty(index_c) || length(index_c) <= 3
            mean_connect(index_x) = mean_liu(phase_y(index_s));        
            else
            mean_connect(index_x) = mean_liu(phase_y(index_c));
            end
        %
        end    
%--------------------------------------------------------------------
%   Start line shift by unwrapping the means    
%--------------------------------------------------------------------
%   unwrap the mean values before global shfit the phase data
    mean_u = mean_connect(1:xr_tmp);
    mean_unwrap = zeros(1,xr_tmp);
    mean_unwrap(1:xr_tmp) = unwrap(mean_u);
%----------------------------------------------------------
    mean_unwrap(1:xr_tmp) = mean_unwrap(1:xr_tmp)...
                                         - round(mean_unwrap(cen_line)/pi_2)*pi_2;
%--------------------------------------------------------------------------
%   shift the phase data
    for index_x = 1:xr_tmp
        index_s = cell_signal{1,index_x};        
        diff_test = mean_unwrap(index_x) - mean_connect(index_x);
        if abs(diff_test) > pi
        phase_tmp(index_s,index_x) = phase_tmp(index_s,index_x)...
                               + pi_2*round(diff_test/(pi_2));                                      
        end
    end
%
    %
out_2D = phase_tmp;
%--------------------------------------------------------------------------

