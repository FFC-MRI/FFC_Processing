% load the dataset using FFCprocessing
% export the FFC acquisition to the workspace
% save in the corresponding analysis folder
% run this script at its location

clear
root = cd;
rawDataFolder = '.\copy of raw data'; % folder that contains the raw data
processedData = '.\test'; % folder where the processed data is saved
cd(rawDataFolder)
folders = dir;
sample = {};
for fi = 3:numel(folders)
    f = folders(fi);
    cd(root)
%     if exist(f.name,'dir')
%         continue
%     end
    
    % to be commented for complete analysis
    if ~isequal(f.name,'TBID_9012 240822')
        continue
    end

    cd(rawDataFolder)
    % loading the data
    try
        cd(f.name)
        load PulseSequenceList.mat
    catch
        cd ..
        continue
    end
    seqs = cellfun(@class,saveList,'UniformOutput',false); % find the sequence names
    indffc = find(cellfun(@(x) strcmp(x,'H9_se_nav_v9'),seqs)); % finds the FFC sequences
    if isempty(indffc)
        continue % skip if no FFC sequence found
    end
    skipflag = 1; % initialise a skip flag
    for indComplete = indffc(end:-1:1) % find the first complete dataset
        kspace = saveList{indComplete}.data; % to be compared with manual despiking below.
        % despiking
%         [kspace,~] = autodespike(kspace,60,6);
        switch f.name % manual deskpiking (VM)
            case 'TBID_8931 040722'
                kspace(37:67,56,1,1,1,2,1) = 0;
            case 'TBID_9012 240822'
                kspace(6:32,23,1,1,1,5,1) = 0;
                kspace(38:90,57,1,1,1,1,3,3) = 0;
            case 'TBID_9122 160123'
                 kspace(7:60,2,1,1,1,5,2) = 0;
                 kspace(2:12,39,1,1,1,2,3) = 0;
                 kspace(33:80,25,1,1,1,4,2) = 0;
            case 'TBID_9112 090123'    
                 kspace(53:90,58,1,1,1,3,2) = 0;
            case 'TBID_8817 130422'
                 kspace(:,27,1,1,1,5,3) = 0;
            case 'TBID_9113 090123'
                 kspace(4:45,45,1,1,1,1,2) = 0;    
            case 'TBID_9094 011222'
                 kspace(69:86,55,1,1,1,5,3) = 0;  
            otherwise
%                  continue % skip processing if the folders dont have spikes
        end
        sze = size(kspace);
        if length(sze)==8 % check that there is only one channel (specific to prostate images)
%             continue
            kspace = kspace(:,:,:,:,:,:,:,1);
            sze = size(kspace);
        end
        if length(sze)>=7 % chec ksize consistency, to avoid processing images that were aborted
            if isequal(sze(6:7),[5,4])||isequal(sze(6:7),[6,4])
                skipflag = 0;
                break
            end
        end
    end
    if skipflag
        continue
    end

    %
    nlines = sze(2);
    nB0 = sze(7);
    nTe = sze(6);
    % estimating the frequency shift
    fids = kspace(:,:,:,:,2,:,:);
    fids(size(fids,1)*10,:) = 0; % oversampling
    spect = fft(fids,[],1);
    [~,ind] = max(abs(fftshift(spect,1)));
    for imchange = 2:nB0*nTe % remove the frequency steps due to changes in image contrast
        ind(1,:,1,1,1,imchange:end) = ind(1,:,1,1,1,imchange:end) + (ind(1,end,1,1,1,imchange-1)-ind(1,1,1,1,1,imchange));
    end
    indref = mean(ind(:)); % reference frequency
    for l = 1:nlines
        for t = 1:nTe
            for b = 1:nB0
                kspace(:,l,1,1,1,t,b) = kspace(:,l,1,1,1,t,b)*exp(-1i*(ind(1,l,1,1,1,t,b)-indref)/10);
            end
        end
    end
    % phase-encoding artefact correction. 
    cd(root)
    if ~exist(f.name,'dir')
        mkdir(root, f.name)
    end
    cd(f.name)
    if isfile('mask.mat')
        load mask.mat
        bkgd = ~mask;
    else
        bkgd = [];
    end
    [imcor,ph,~,bkgd,jj] = iterative_images_correction_v7(kspace(:,:,1,1,1,:,:),0.15,20,0.15,bkgd);  
    % saving the results
    save initial_image.mat imcor bkgd ph kspace
    if ~isfile('mask.mat')
        mask = sum(~bkgd,3)>(size(bkgd,3)/8);
        save mask.mat mask
    end
    
    %----------------------------------------
    % modification VM 02/05/23
    imcor = squeeze(imcor)/max(imcor(:));
    N_times = size(imcor,3);
    N_fields = size(imcor,4);
    adjustment_factor = 1.2;
    bkgd = ~mask;  % assess the background region
    imcordenoised = zeros(size(imcor));

     for i = 1:N_fields
         for  j = 1:N_times
             img = imcor(:,:,j,i);
             noise = adjustment_factor*std(img(bkgd(:)));
             imcordenoised(:,:,j,i) = BM3D(real(img),noise) + 1i*BM3D(imag(img),noise);
        end
     end
    save denoised_image.mat imcordenoised
    %---------------------------------------
    
    % get the T1 maps
    B = [];
    T = [];
    for i = 1:length(saveList{indComplete}.waveformProfile.waveList)
        B(i) = saveList{indComplete}.waveformProfile.waveList{i}.Bevo;
        T(i) = saveList{indComplete}.waveformProfile.waveList{i}.Tevo;
    end
    [B,elnum] = unique(B);
    [~,experimentOrder] = sort(elnum);
    B = B(experimentOrder);
    
    % make sure the dispersion profile is provided correctly
    [B,order] = sort(B,'descend');
    B = B';
    T = reshape(T,sze(end-1),sze(end))';
    T = T(order,:);
    imcordenoised = imcordenoised(:,:,:,order);
    B = B*1000; % convert into mT
    T = T*1000; % convert into ms
    outfit = run_fit(abs(imcordenoised),T,B,200,mask,0,0);
%     outfit.R1highres = imresize(outfit.R1,2,"bilinear"); % save a high-resolution version
    save fitresult.mat outfit
    
    % add the maps to the list of data
    outfit.patient = f.name;
    sample{end+1} = outfit;
end
cd(root)

bkgd = squeeze(bkgd(:,:,1));
% figure;imagesc(abs(imcordenoised(:,:,1,4).*~bkgd))
% axis off;colormap(gray)
% figure;imagesc(abs(imcor(:,:,5,1).*~bkgd))
% axis off;colormap(gray)
% figure;imagesc(outfit.R1highres(:,:,4));
% axis off
figure;imagesc(outfit.R1(:,:,1).*~bkgd);
axis off;
% figure;imagesc(outfit.rsquare);
% axis off