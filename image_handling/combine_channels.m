function [combined_images] = combine_channels(images,opts,obj)
%multicoil reconstruction
%   Detailed explanation goes here

tic
tempdim = size(images);

images = reshape(images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints,obj.n_receivers]); % old debug code, enforce a specific dimensionality (with singletons if needed)

images_temp = reshape(images,tempdim(1),tempdim(2),[],obj.n_receivers); %reduce everything to 4 dimensions (X,Y, everything else, number of receivers).


switch opts
    case 1
        n_channels = size(images,6);
        if n_channels >1
            for s=1:size(images_temp,3)
                
                psi = [n_channels,n_channels];
                dims = size(images_temp);
                channel_data = reshape(images_temp(:,:,s,:),[],n_channels);          %make everything into a vector so we don't have to worry about every possible combination of slices, echoes etc
                eta1 = (images_temp(:,:,s,:));
                eta = squeeze(eta1(1,:,1,:))'; %Assume the first line of k-space is mostly noise. Almost universally the case at 0.2T.
                psi = (1/(length(eta)-1))*(eta*eta');
                try
                    L = chol(psi,'lower');
                    L_inv =inv(L);
                catch
                    L_inv =1; %This usually indicates something bad has happened but prevents complete failure
                end
                %                             L_inv =1;
                
                data_scaled = ((L_inv)*(permute(channel_data,[2 1]))); %scale the data depending on the noise figure from each element.
                
                channel_data = reshape(data_scaled',tempdim(1),tempdim(2),[]);
                combined_images(:,:,s) = sqrt(sum(channel_data.*conj(channel_data),3)); %equivilent to rssq() I think. Tutorial is explicit about using the conjugate for reasons I don't understand
                
            end
            
            combined_images = reshape(combined_images,[tempdim(1),tempdim(2),obj.slices,obj.n_timepoints,obj.n_fieldpoints]); %put everything back to its original dimensionality.
            
        else
            combined_images = mean(images,6); %If something is wrong and we reach this point, just average over the multicoil dimension.
        end
    case 2
        combined_images = mean(images,6); %User has selected to average over the multicoil data. Not usually a good idea but useful in particular conditions.
    otherwise
        combined_images = images(:,:,:,:,:,opts-2); %displays the individual channel data for debugging purposes.
        
        
end

    function S = adaptive_est_sens(data) %adapted from Oxford SENSE tutorial based on Walsh et al for testing purposes.
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



