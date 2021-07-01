function [combined_images] = combine_channels(images,noise,opts,obj)
%multicoil reconstruction
%   Detailed explanation goes here
%
% noise = reshape(noise,numel(noise)/2, 2);
% noise = permute(noise,[2 1]);
% M = size(noise,2);
% Rn = (1/(M-1))*(noise*noise');
% Bn = (mean(abs(noise(:)))^2)./abs(noise(length(noise)./2))^2;
% Rnscaled = (Rn);

tempdim = size(images);

images = reshape(images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);

switch opts
    case 1
        n_channels = size(images,6);
        if n_channels >1
            for s=1:size(images,3)
                images_temp = images(:,:,s,:,:,:);
                eta = squeeze(noise(:,s,:));
                psi = [n_channels,n_channels];
                dims = size(images_temp);
                channel_data = (reshape(images_temp,[],n_channels));
                psi = (1/(length(squeeze(noise(:,s,:)))-1))*(eta*eta');
                try
                    L = chol(psi,'lower');
                    L_inv =inv(L);
                catch
                    L_inv =1;
                end
                %    L_inv =1;
                data_scaled = ((L_inv)*permute(channel_data,[2 1]));
                channel_data = reshape(data_scaled',dims);
                combined_images(:,:,s,:,:) = sqrt(sum(channel_data.*conj(channel_data),6));
%                 combined_images = rssq(images,6); %basic S
            end
        else
            combined_images = mean(images,6);
        end
    case 2
        combined_images = mean(images,6);
    otherwise
        combined_images = images(:,:,:,:,:,opts-2);
end

% for n=1:size(images,6)
% coil_images(:,:,:,n) = coil_images(:,:,:,n).*(noise(1)./noise(n));
% end
%
% combined_images = sqrt(sum(coil_images.*conj(coil_images),4));