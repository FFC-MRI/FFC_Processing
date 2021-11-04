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
tic
tempdim = size(images);

images = reshape(images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]);

images_temp = reshape(images,tempdim(1),tempdim(2),[],obj.n_receivers);


switch opts
    case 1
        n_channels = size(images,6);
        if n_channels >1
            for s=1:size(images_temp,3)
                
                eta = squeeze(noise(:,1,:));
                psi = [n_channels,n_channels];
                dims = size(images_temp);
                channel_data = reshape(images_temp(:,:,s,:),[],n_channels);
                
                %                  channel_data(:,2) = channel_data(:,2).*exp(1i*-pi);
                
                
                
                eta1 = (images_temp(:,:,s,:));
                eta = squeeze(eta1(1,:,1,:))'; %Assume the first line of k-space is mostly noise
                psi = (1/(length(eta)-1))*(eta*eta');
                try
                    L = chol(psi,'lower');
                    L_inv =inv(L);
                catch
                    L_inv =1; %This usually indicates something bad has happened but prevents complete failure
                end
                %                             L_inv =1;
                
                data_scaled = ((L_inv)*(permute(channel_data,[2 1])));
                
                channel_data = reshape(data_scaled',tempdim(1),tempdim(2),[]);
                combined_images(:,:,s) = sqrt(sum(channel_data.*conj(channel_data),3));
                
                %                 combined_images(:,:,s) = sqrt(sum(abs(channel_data).^2,3));
                
                %                  kernel = ones(27)*0.00000001;
                
                
%                 sensitivities = squeeze(adaptive_est_sens(reshape(channel_data,[tempdim(1),tempdim(2) 1 n_channels])));
                
                
                %
                %
                %
                img_combined = combined_images(:,:,s);
                % %                 thresh = 0.05*max(abs(img_combined(:)));
                % % mask = abs(img_combined) > thresh;
                % % sensitivities = sensitivities.*mask;
                %
%                 img_combined_opt = sum(channel_data.*conj(sensitivities),3)./sum(sensitivities.*conj(sensitivities),3);
                %            [img_combined_opt,smap] = ir_mri_coil_combine(channel_data);
                % %
                % %
                % %
                % %
%                 combined_images(:,:,s) = img_combined_opt;
                %                   combined_images(:,:,s) = rssq(channel_data,3); %basic S
            end
            
            combined_images = reshape(combined_images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints]);
            
        else
            combined_images = mean(images,6);
        end
    case 2
        combined_images = mean(images,6);
    otherwise
        combined_images = images(:,:,:,:,:,opts-2);
        
        
end

    function S = adaptive_est_sens(data) %adapted from Oxford SENSE tutorial based on Walsh et al. Seems to work quite well for low field data
        [Nx,Ny,Nz,Nc] = size(data);
        S = zeros(Nx,Ny,Nz,Nc);
        M = zeros(Nx,Ny,Nz);
        w = 5;
        
        for i = 1:Nx
            ii = max(i-w,1):min(i+w,Nx);
            for j = 1:Ny
                jj = max(j-w,1):min(j+w,Ny);
                for k = 1:Nz
                    kk = max(k-w,1):min(k+w,Nz);
                    kernel = reshape(data(ii,jj,kk,:),[],Nc);
                    [V,D] = eigs(conj(kernel'*kernel),1);
                    S(i,j,k,:) = V*exp(-1j*angle(V(1)));
                    M(i,j,k) = sqrt(D);
                end
            end
        end
        
        %      S = S.*(M>0.1*max(abs(M(:))));
    end
end




% for n=1:size(images,6)
% coil_images(:,:,:,n) = coil_images(:,:,:,n).*(noise(1)./noise(n));
% end
%
% combined_images = sqrt(sum(coil_images.*conj(coil_images),4));