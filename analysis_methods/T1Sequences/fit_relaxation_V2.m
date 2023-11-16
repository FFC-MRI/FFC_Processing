classdef fit_relaxation_V2 < R1Template
%Give the Name that will be compared to in the GUI and Checker
    properties
    processName = 'Fit Relaxation V2';
    end


%Example Code
    methods
        function [R1out,others,tres,sse,err,fitres] = FitRelaxation(imageSmoothed,t,B0,B0_pol,rois)
            %% initialisation
            [number_fields, number_times] = size(t);
            t = t./1000; % convert to sec
            sze = size(imageSmoothed);
            R1out = zeros(sze(1)*sze(2),number_fields);
            others = zeros(sze(1)*sze(2),3);
            tres = zeros(sze(1)*sze(2),1);
            sse = zeros(sze(1)*sze(2),1);
            err = zeros(sze(1)*sze(2),number_fields);
            roiInd = find(rois(:)>0); % position of the voxels to be analysed in the image

            %% data vectorisation: needed for efficient parfor loops.
            imageVextorised = reshape(squeeze(imageSmoothed),sze(1)*sze(2),number_times,number_fields);
            imageVextorised = imageVextorised(roiInd,:,:);
            R1outVect = R1out(roiInd,:);
            othersVect = others(roiInd,:); 
            tresVect = tres(roiInd); 
            sseVect = sse(roiInd);
            errVect = err(roiInd,:);

            %% setting up surface fit 
            R1name = cell(numel(B0),1);
            h = cell(numel(B0),1);
            R1list = '';
            % generate the list of relaxation rate variables
            for nb = 1:numel(B0)
                R1name{nb} = ['R1_' num2str(nb)];
                R1list = [R1list R1name{nb} ', ']; %#ok<*AGROW>
                R1Upper(nb) = 150;
                R1Lower(nb) = 0.3;
                R1Start(nb) = 1.2*B0(nb)^-0.3;
            end
            % generate the model for each field
            modelStr = [];
            for nb = 1:numel(B0)
                h{nb} = singleexphandle_ns(R1name{nb},B0(nb));
                modelStr = [modelStr 'h{' num2str(nb) '}(' R1name{nb} ',a_dw, b_dw, chi,t,B0) + '];
            end
            % assemble the models into one function
            eval(['fun = @(' R1list 'a_dw, b_dw, chi,t,B0) ' modelStr(1:end-2) ';'])
            ftype = fittype(fun, 'independent', {'t', 'B0'});

            %% setting up fit parameters, x and y variables (these are the same for all voxels)
            x = t(:); 
            y = repmat(B0,1,number_times); 
            y = y(:);

            %% fitting loop
            parfor ix = 1:length(roiInd) 
                % preparing the loop-dependant parameters and variables
                foption = fitoptions(ftype);
                foption.Robust = 'Bisquare';
                foption.MaxIter = 10000;
                foption.MaxFunEvals = 10000;
                z = squeeze(imageVextorised(ix,:,:))';
                z = z(:); 
                chiGuess = mean(z(B0_pol>0,end)./B0_pol); % estimation for the magnetic susceptibility
                foption.Upper = [R1Upper 10*B0_pol 10*B0_pol chiGuess*100];  % R1_1,R1_2,R1_3,R1_4,Bpol*alpha, Bread*beta, chi
                foption.Lower = [R1Lower 0         0         0];
                foption.StartPoint = [R1Start B0_pol*0.5 B0_pol*0.1 chiGuess];
                % perform the fit and store the results
                [fitobject,gof] = fit([x,y],z,ftype,foption);
                bestFit = coeffvalues(fitobject);
                R1outVect(ix,:) = bestFit(1:numel(B0));
                othersVect(ix,:) = bestFit((numel(B0)+1):end);
                tresVect(ix) = gof.rsquare; 
                sseVect(ix) = gof.sse;
                % extracting the confidence intervals:
                ci = confint(fitobject);
                errVect(ix,:) = diff(ci(:,1:numel(B0)),1,1);
                fitres{ix} = fitobject;
            end

            %% formats the vectorised outputs back into images
            R1out(roiInd,:) = R1outVect;
            others(roiInd,:) = othersVect; 
            tres(roiInd) = tresVect; 
            sse(roiInd) = sseVect;
            err(roiInd,:) = errVect;

            R1out = reshape(R1out,sze(1),sze(2),number_fields);
            others = reshape(others,sze(1),sze(2),size(othersVect,2));
            tres = reshape(tres,sze(1),sze(2));
            sse = reshape(sse,sze(1),sze(2));
            err = reshape(err,sze(1),sze(2),number_fields);

        end
        
        function h = singleexphandle_ns(R1name,Bevo)
        % TODO: eval is bad, it should be changed for feval or similar. 
        if Bevo < 0.1  % prepolarised
            eval(['h = @(' R1name ', a_dw, b_dw, chi, t, B0) model(' R1name ', a_dw, b_dw, chi, t, B0);'])
        else  % non-polarised
            eval(['h = @(' R1name ', a_dw, b_dw, chi, t, B0) model(' R1name ', 0, b_dw, chi, t, B0);'])
        end

            function Mz = model(R1, a_dw, b_dw, chi, t, B0)
                Mz = (a_dw.*exp(-t*R1) + ...        % contribution from the polarisation
                      B0.*(1-exp(-t*R1)) + b_dw) ...     % contribution from the evolution field
                     .*(B0==Bevo)*chi;      % restriction to the field of interest and scaling

    %             Mz = ((chi*Bpol*alpha).*exp(-t*R1) + ...        % contribution from the polarisation
    %                   (chi*Bevo - MrampUp).*(1-exp(-t*R1))) ...     % contribution from the evolution field
    %                  .*(B0==Bevo)*chi;      % restriction to the field of interest and scaling

    %             Mz = ((eta*chi*Bpol + a*B0    + b).*exp(-t*R1) + ...        % contribution from the polarisation
    %                   (chi*B0       + a*Bread + b).*(1-exp(-t*R1))) ...     % contribution from the evolution field
    %                  .*(B0==Bevo);      % restriction to the field of interest
    %           a_dw = a
    %           b_dw = eta*chi*Bpol + b
    %           a_up = chi
    %           b_up = a*Bread + b
                end
        end
        
    end
end

