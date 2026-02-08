function [filtered_images] = ffc_mri_filter(images, filter_type, kernel)
%FFC_MRI_FILTER  Drop-in denoising pipeline for magnitude stacks (low-field friendly)
%
% images: [X Y ...] (any additional dims flattened)
% filter_type:
%   'LLR+Wavelet'        (recommended default)
%   'LLR'               (faster)
%   'Wavelet'           (fast baseline)
%   'LLR+Wavelet+TGV'    (strongest, slowest)
%   others: passthrough magnitude
%
% kernel: kept for compatibility; unused here except your other cases.

dim = size(images);
temp = abs(reshape(images, dim(1), dim(2), []));   % [X Y N]

switch filter_type
    case 'LLR+Wavelet'
        temp = llr_denoise_stack(temp, struct('patch',8,'step',3,'niter',2,'tau',2.2,'sigma',[]));
        temp = wavelet_per_frame(temp, 3); % level 3 for 128x128

    case 'LLR'
        temp = llr_denoise_stack(temp, struct('patch',8,'step',3,'niter',2,'tau',2.2,'sigma',[]));

    case 'Wavelet'
        temp = wavelet_per_frame(temp, 3);

    case 'Total Variation'
        temp = llr_denoise_stack(temp, struct('patch',8,'step',3,'niter',2,'tau',1.2,'sigma',[]));
        temp = wavelet_per_frame(temp, 3);
        % very light TGV (optional but effective on residual structured junk)
        optsT = struct('data','L1','normalize',true);
        tmp2 = zeros(size(temp), 'like', temp);
        parfor n = 1:size(temp,3)
            tmp2(:,:,n) = TGV2denoise(temp(:,:,n), 0.04, 0.015, 25, 150, optsT);
        end
        temp = tmp2;

    otherwise
        % passthrough magnitude
end

filtered_images = reshape(temp, dim);
end
