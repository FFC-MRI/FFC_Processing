function [stackCorrected, phase, bkgd, iterationUsed, exitflag] = iterative_images_correction_v8(stackFft, thresh_bkgd, max_iterations, thresh_phase, known_bkgd, regFlag)

% iterative_images_correction_v8
%
% [stackCorrected,phase, bkgd, iterationUsed, exitflag] = 
% iterative_images_correction_v8(stackFft, thresh_bkgd, max_iterations,
% thresh_phase, known_bkgd, regFlag) 
%
% Inputs:
% stackFft: multi-dimensional image stack, k-space. The only requirement is
% to set dimension 1 as frequency encode, 2 as phase encode and last one as
% receiver number. Make sure that the k-space is centred (use FFTSHIFT if
% not).
% thresh_bkgd: value used to threshold the background in the reconstructed
% image. 
% max_iterations: maximum number of iterations to use.
% thresh_phase: threshold for the phase correction, to decide that we
% reached convergence.
% known_bkgd: 2D mask situated in a region of the image known to be empty.
% regFlag: flag for image registration (1 to register imaged, 0 otherwise)
%
% Typical values: thresh_bkgd = 0.15;max_iterations = 20;thresh_phase =
% 1e-6; known_bkgd =[];regFlag = 0;    
% 
% Outputs:
% stackCorrected: stack of corrected images, with the same dimensionality
% as the input stack.
% phase: estimation of the phase error for each image.
% bkgd: estimation of the background.
% iterationUsed: number of iterations required for convergence.
% exitflag: 1 if convergence was reached, 0 otherwise.
%
% Author: Lionel Broche, University of Aberdeen, March 2017
% License: GPL v3 (https://www.gnu.org/licenses/gpl.html)

% starts by making sure that both the image and the Fourier space are
% centric
stackFft = RecentreImage(stackFft);
stackCorrected = zeros(size(stackFft));

if ndims(stackFft)>3
    % case of a multi-dimentional matrix. Multi-slice images must be
    % treated separately as the background differs between slices
    if size(stackFft,4)>1 % if there are more than one slice
        for slice = 1:size(stackFft,4)
            % isolate the correct background
            if ~isempty(known_bkgd)
                if ndims(known_bkgd)==3
                    known_bkgd_slice = squeeze(known_bkgd(:,:,slice));
                else
                    known_bkgd_slice = known_bkgd;
                end
            else
                known_bkgd_slice = known_bkgd;
            end
            [stackCorrected(:,:,:,slice,:,:,:,:,:), phase(:,:,slice,:,:,:,:,:), bkgd(:,:,slice), iterationUsed(slice), exitflag] = ...
                iterative_images_correction_v8(stackFft(:,:,:,slice,:,:,:,:,:), thresh_bkgd, max_iterations, thresh_phase, known_bkgd_slice, regFlag);
        end
        return
    end
    % Start by finding the image with
    % best SNR to get a good estimation of the background for the others.
    % this 'if' section splits the multidimensional matrix into 2 or 3D
    % elements that are processed recursively.
    sze = [1 1 1 1 1 1 1 1 1];
    sze(1:ndims(stackFft)) = size(stackFft); % dimensions are: signal, line, views2, slice, echoes, average, time, field, channel
    stackFft = reshape(stackFft,sze(1),sze(2),prod(sze(3:end-1)),sze(end));
    stackCorrected = zeros(size(stackFft));
    phase = zeros(sze(2),size(stackFft,3));
    meanSignal = squeeze(mean(median(median(abs(stackFft),1),2),4));
    [~,processOrder] = sort(meanSignal,'descend');
    stackFft = stackFft(:,:,processOrder,:);
    if regFlag
        % registration over the frequency-encode direction only. Quite
        % bulky, this may be optimised a lot.
        % check register_images, modify for 1D registration
        upsamp = 10; % sub-voxel resolution factor
        stack = fft(fft(stackFft,[],1),[],2); % using the real space
        regSig = squeeze(sum(sum(abs(stack),2),4)); % we only consider the average signal over the freq-encode direction, to avoid smearing
        regSigFft = fft(regSig,[],1); % we find the sub-voxel shift by comparing how the object move in the freq-encode direction
        regSigFftPad = zeros(upsamp*size(regSigFft,1),size(regSigFft,2));
        regSigFftPad((1:size(regSigFft,1))+size(regSigFft,1),1:size(regSigFft,2)) = fftshift(regSigFft,1);
        regSigPad = abs(ifft(regSigFftPad));
        phShift = zeros(size(stackFft));
        for i = 1:size(regSigPad,2)
            regSigPad(:,i) = regSigPad(:,i)/max(regSigPad(:,i)); % normalisation
            regSigPad(:,i) = regSigPad(:,i) > 0.25; % masking, to avoid prolems due to T1 contrast
        end
        for i = 2:size(regSigPad,2) % shifting the images
            c = xcorr(regSigPad(:,1),regSigPad(:,i));
            voxel = ((1:length(c)) - length(c)/2)/upsamp;
            [~,ind] = max(c);
            dv(i) = voxel(ind);
            ph = linspace(-pi,pi,size(stackFft,1))'*dv(i);
            ph = ph*ones(1,size(stackFft,2));
            for rec = 1:size(phShift,4)
                stackFft(:,:,i,rec) = stackFft(:,:,i,rec).*exp(1i*ph);
            end

        end
    end
    % correct the image with highest SNR first, using all channels
    [stackCorrected(:,:,1,:),phase(:,1),bkgd,iterationUsed,exitflag(1)] = ...
        iterative_images_correction_v8(squeeze(stackFft(:,:,1,:)), thresh_bkgd, max_iterations, thresh_phase, known_bkgd,0);
    % then use the background estimation to correct the other images
    parfor i = 2:size(stackFft,3)
        [stackCorrected(:,:,i,:),phase(:,i),bkgdStack(:,:,i),~,exitflag(i)] = ...
            iterative_images_correction_v8(squeeze(stackFft(:,:,i,:)),thresh_bkgd,1,thresh_phase,bkgd,0);
%             iterative_images_correction_v8(squeeze(stackFft(:,:,i,:)),thresh_bkgd,1,thresh_phase,bkgd0,1); % only use one iteration, do not optimise the background (TO FIX)
    end
    % there is a bit or re-ordering to be done here
    [~,reverseOrder] = sort(processOrder,'ascend');
    stackCorrected = stackCorrected(:,:,reverseOrder,:); % re-order for the pre-correction sorting
    phase = phase(:,reverseOrder);
    % then re-shape into the original size
    stackCorrected = reshape(stackCorrected,sze);
    phase = reshape(phase,sze(2:end-1));
    bkgd = bkgd;
    return
end

% if the program gets here, then we have either a 2D or 3D image. 3D images
% are normal 2D images from several channels. To make things easier, we
% reshape 2D images into 3D.
if ismatrix(stackFft)
    stackFft = reshape(stackFft,size(stackFft,1),size(stackFft,2),1);
end

% reconstructs the image for each channel:
stack = ifft(ifft(stackFft,[],1),[],2);
bkgd = false(size(stack));

% generates the mask
if isempty(known_bkgd)
    known_bkgd = false(size(stackFft,1),size(stackFft,2));
%     level = graythresh(mean(abs(stack),3));
%     bkgd = im2bw(mean(abs(stack),3),level); % fails, the image shifts a lot

    maxValue = max(max(abs(stack),1),2);
    for channel = 1:size(stack,3)
        bkgd(:,:,channel) = abs(stack(:,:,channel))<(thresh_bkgd*maxValue(channel)); % Method using the magnitue data
    end
%     for channel = 1:size(stack,3)
%         bkgd(:,:,channel) =
%         stdfilt(angle(squeeze(mean(stack(:,:,channel),3))))>1;  % method using the phase data. Not robust...
%     end
else 
    bkgd = known_bkgd==1;
end

% initialisation:
ph0 = zeros(1,size(stackFft,2));       % initial guess for the phase. Setting it to 0 allows to get the correction increment at each step, which should tend to zero.
stackCorrectedFft = stackFft;           % initial data is the damaged data
phase = zeros(1,size(stackFft,2));        % list of phase estimates

% iterative correction loop:
for iterationUsed = 1:max_iterations
    [dph,stackCorrectedFft,stackCorrected,fval,exitflag] = optim_image_v8(stackCorrectedFft,bkgd,ph0,max_iterations == 1);    % optimises the phases
    phase = phase + dph;              % cumulate the phase corrections
    
    if (std(dph)<thresh_phase)||(iterationUsed>=max_iterations)    % convergence is reached when the amplitude of the phase correction falls below the threshold
        break
    end
    % re-defines the background, remove the first line of the k-space in case of DC artefact
%     maxValue = max(max(abs(stackCorrected),1),2);
%     for channel = 1:size(stackCorrected,3)
%         bkgd(:,:,channel) = (known_bkgd(:,:,channel)==1)|abs(stackCorrected(:,:,channel))<(thresh_bkgd*maxValue(channel));
%     end
    
    for channel = 1:size(stackCorrected,3)
%         bkgd(:,:,channel) = stdfilt(angle(squeeze(mean(stackCorrected(:,:,channel),3))))>1;
        bkgd(:,:,channel) = (known_bkgd==1) | (abs(stackCorrected(:,:,channel))<(thresh_bkgd*max(max(abs(stackCorrected(2:end,2:end,channel))))));
    end    
    bkgd = squeeze(sum(bkgd,3))==channel;

end

% 
% % avoid check board pattern problem in phase space
% stackCorrectedFft = fftshift(fftshift(stackCorrectedFft,1),2);
% stackCorrected = fftshift(fftshift(ifft(ifft(stackCorrectedFft,[],1),[],2),1),2);
