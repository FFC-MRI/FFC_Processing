function [phaseCorrected,stackCorrectedFft,stackCorrected,fval,exitflag] = optim_image_v8(stackFft,bkgd,phaseInitial,preciseFlag)
% corrects for phase encode artifacts due to random phase fluctuations
% between k-space lines
% [phaseCorrected,stackCorrectedFft,stackCorrected,fval,exitflag] = optim_image_v8(stackFft,bkgd,phaseInitial,preciseFlag)
% inputs: 
% stackFft: 3D image, third dimension is the recording channel number.
% the phase encode direction should be vertical in the image (along 
% columns)
% bkgd: 3D image, estimation of the background for all images.
% phaseInitial: initial estimation of the phase error.
% preciseFlag: 0 for fast calculation, likely to fail to converge but helps
% correcting the background; 1 for long calculation that are likely to
% converge (when the background is known)
% outputs:
% phaseCorrected: array, estimation of the phase error along each column
% stackCorrectedFft: 3D image, corrected k-space
% stackCorrected: 3D image, corrected images in real space (complex)
% fval: min value returned by fminsearch. Really useless.
% exitflag: exit status provided by fminsearch.
%
% Lionel Broche, 08/11/2016

% the initial stack of image is of dimension 3, with the 3rd dimension
% being the number of acquisition channels
[l,c,channelNum] = size(stackFft);
% stackFftNorm = stackFft/max(stackFft(:)); % normalise the image for standardised calculation

% optimisation parameters
if preciseFlag
    options = optimoptions('fmincon','Algorithm','sqp',...
                             'MaxFunEvals',5e4,...
                             'MaxIter',100000,...
                             'Display','off',...
                             'TolX',1e-7,...
                             'UseParallel',true);
     % initialising variables
    Aeq = zeros(1,c);   % fmincon satisfies the condition Aeq*x = beq
    beq = zeros(1,1); 
    beq(1) = 0; Aeq(1,1:c) = 1; % The average phase correction is set to 0, as this parameter is free
else
    options = optimoptions('fmincon','Algorithm','sqp',...
                             'MaxFunEvals',1e3,...
                             'MaxIter',5000,...
                             'Display','off',...
                             'TolX',1e-2,...
                             'UseParallel',true);
     % initialising variables
    Aeq = zeros(2,c);   % fmincon satisfies the condition Aeq*x = beq
    beq = zeros(1,2); 
    beq(1) = 0; Aeq(1,1:c) = 1; % The average phase correction is set to 0, as this parameter is free
    beq(2) = 0; Aeq(2,1:c) = 1:c; % Set the first moment to zero to avoid image shifting
end
                   
[phaseCorrected,fval,exitflag] = fmincon(@(x) minabsim(x),phaseInitial',[],[],Aeq,beq,-10*pi*ones(1,c), 10*pi*ones(1,c),[],options); % corrects for random phase on all lines

% corrects the image using the latest estimation of the phase error
phaseCorrected = linear_correction(phaseCorrected)';
correctionMatrix = repelem(ones(l,1)*exp(-1i*phaseCorrected),1,1,channelNum);
stackCorrectedFft = stackFft.*correctionMatrix;
stackCorrected = ifft(ifft(stackCorrectedFft,[],1),[],2);

    function s = minabsim(phi)
        correctionMatrix = repelem(ones(l,1)*exp(-1i*phi'),1,1,channelNum); % repeat the phase correction for each channel
        % as fft2 does not work well on 3D stacks of 2D images, we use fft twice instead
        Ac = stackFft.*correctionMatrix;   % apply the correction in the k-space
        Ic = ifft(ifft(Ac,[],1),[],2);       % get the corrected image
        Ic = abs(Ic);              % only consider the absolute value
        s = sum(Ic(bkgd(:)));               % assess the amount of noise in the background
    end

end
