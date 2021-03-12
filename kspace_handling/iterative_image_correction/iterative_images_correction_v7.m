function [I1,ph,A1,bkgd,jj] = iterative_images_correction_v7(A,thresh_bkgd,max_iterations,thresh_phase,known_bkgd)

% [I1,ph,A1,bkgd,jj] = iterative_images_correction_v7(A,thresh_bkgd,max_iterations,thresh_phase,known_bkgd)
%
% Typical values: thresh_bkgd = 0.15;max_iterations = 20;thresh_phase = 0.15; known_bkgd =[];
%
%

I1 = [];
ph = [];
A1 = [];
bkgd = [];
jj = [];

if ndims(A)>2
    sze = size(A);
    bkgd = [];
    A = reshape(A,size(A,1),size(A,2),[]);


    parfor n = 1:size(A,3)
        [I1(:,:,n),ph(:,:,n),~,~,~] = iterative_images_correction_v7(A(:,:,n),thresh_bkgd,max_iterations,thresh_phase,known_bkgd);
    end
     I1 = reshape(I1,sze);
     return
end


% reconstructs the image:
a = ifft2c(A);

% generates the mask
if nargin<=4
    known_bkgd = zeros(size(A))==1;
    bkgd = abs(a)<(thresh_bkgd*max(abs(a(:))));
elseif isempty(known_bkgd)
    known_bkgd = zeros(size(A))==1;
    bkgd = abs(a)<(thresh_bkgd*max(abs(a(:))));
else
    bkgd = known_bkgd==1;
end

% initialisation:
ph0 = zeros(1,size(A,2));       % initial guess for the phase. Setting it to 0 allows to get the correction increment at each step, which should tend to zero.
A1 = A;                         % A1 is the corrected k-space. This only allocates enough memory
I1 = a;                         % same for the corrected image
ph = zeros(1,size(A,2));        % list of phase estimates
bkgd_score = 0;

% iterative correction loop:
for jj = 1:max_iterations
    [dph,A1,I1,fval,exitflag] = optim_image_v7(A1,bkgd,ph0);     % optimises the phases
    ph = ph + dph;              % cumulate the phase corrections
    %     % analysis of the histogram to extract an estimation of the background
    %     [v,b] = hist(real(I1(:)),1000);
    %     sigma = std(real(I1(bkgd(:))));       % estimation of the noise level
    
    if std(dph)<thresh_phase    % convergence is reached when the amplitude of the phase correction falls below the threshold
        break
    end
    bkgd = (known_bkgd==1) | (abs(I1)<(thresh_bkgd*max(max(abs(I1(2:end,2:end))))));        % re-defines the background, remove the first line of the k-space in case of DC artefact
    %     new_score = sum(bkgd);
    %     if abs(new_score-bkgd_score)<thresh_phase
    %         break
    %     else
    %         bkgd_score = new_score;
    %     end
end


