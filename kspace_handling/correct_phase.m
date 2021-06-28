function [correctedkspace] = correct_phase(kspace,backgroundtest,n_receivers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

dims = size(kspace); %we can reshape back to this dimensionality later

    workingkspace = reshape(kspace(:,:,:,:,:,n_receivers),dims(1),dims(2),[]); %only do the phase correction for one coil and duplicate the results to others.
    workingdim = size(workingkspace);
    
    
    if backgroundtest ==1
        hh = figure;
        imagesc(abs(fft2c(kspace(:,:,1,1,1)))); axis off square; colormap('gray');
        background = roipoly;
        close(hh);
    else
        background = [];
    end
    [I1,ph] = iterative_images_correction_v7(workingkspace(:,:,:,:,:,1),0.20,10,0.00015,background);
    
    phase_correction = reshape(repmat(ph,dims(1),1,1),size(kspace(:,:,:,:,:,1)));
    
    for n=1:n_receivers
        correctedkspace(:,:,:,:,:,n) = kspace(:,:,:,:,:,n).*exp(-1i*phase_correction);
    end
%     
% correctedkspace = kspace; %bypass for debugging

