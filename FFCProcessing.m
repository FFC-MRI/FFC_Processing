classdef FFCProcessing < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        QualityAssurancePanel          matlab.ui.container.Panel
        AllButton                      matlab.ui.control.Button
        THKButton                      matlab.ui.control.Button
        SNRButton                      matlab.ui.control.Button
        UniformityButton               matlab.ui.control.Button
        GeometricAccButton             matlab.ui.control.Button
        MultichannelDropDown           matlab.ui.control.ListBox
        ExportOptionsPanel             matlab.ui.container.Panel
        CreateFigurePanel              matlab.ui.container.Panel
        ImageButton                    matlab.ui.control.Button
        kspaceButton                   matlab.ui.control.Button
        ExportscanparameterstoworkspaceButton  matlab.ui.control.Button
        ExportstackedimagesButton      matlab.ui.control.Button
        ExportallimagesButton          matlab.ui.control.Button
        ExportscantoworkspaceButton    matlab.ui.control.Button
        Panel_2                        matlab.ui.container.Panel
        BackgroundselectionButton      matlab.ui.control.StateButton
        ReprocessingButton             matlab.ui.control.StateButton
        LoadDataSetButton              matlab.ui.control.Button
        RelaxometryPanel               matlab.ui.container.Panel
        ColourAxisLowerLimitEditField  matlab.ui.control.NumericEditField
        ColourAxisLowerLimitEditFieldLabel  matlab.ui.control.Label
        ColourAxisUpperLimitEditField  matlab.ui.control.NumericEditField
        ColourAxisUpperLimitEditFieldLabel  matlab.ui.control.Label
        DisplayingasT1Button           matlab.ui.control.StateButton
        LogarithmicMapsButton          matlab.ui.control.StateButton
        CheckfitsButton                matlab.ui.control.StateButton
        GenerateROIDispersionCurveButton  matlab.ui.control.Button
        GenerateRelaxationMapsButton   matlab.ui.control.Button
        ScanListPanel                  matlab.ui.container.Panel
        ScanList                       matlab.ui.control.ListBox
        DisplayOptionsPanel            matlab.ui.container.Panel
        IndividualNormalisationButton  matlab.ui.control.StateButton
        FlipVerticallyButton           matlab.ui.control.Button
        FlipHorizontallyButton         matlab.ui.control.Button
        RotateImageButton              matlab.ui.control.Button
        WindowWidthSlider              matlab.ui.control.Slider
        WindowWidthSliderLabel         matlab.ui.control.Label
        WindowLevelSlider              matlab.ui.control.Slider
        WindowLevelSliderLabel         matlab.ui.control.Label
        CurrentImagePanel              matlab.ui.container.Panel
        RegisterImagesButton           matlab.ui.control.Button
        ImageTypeDropDown              matlab.ui.control.DropDown
        ImageTypeDropDownLabel         matlab.ui.control.Label
        SliceDropDown                  matlab.ui.control.DropDown
        SliceDropDownLabel             matlab.ui.control.Label
        EvolutionTimeDropDown          matlab.ui.control.DropDown
        EvolutionTimeDropDownLabel     matlab.ui.control.Label
        EvolutionFieldDropDown         matlab.ui.control.DropDown
        EvolutionFieldDropDownLabel    matlab.ui.control.Label
        panel_recon                    matlab.ui.container.Panel
        FFTDropDown_3D                 matlab.ui.control.DropDown
        Sizeof3DFFTLabel               matlab.ui.control.Label
        FilterSettingsButton           matlab.ui.control.Button
        BiasFieldDropDown              matlab.ui.control.DropDown
        BiasFieldDropDownLabel         matlab.ui.control.Label
        CombinationDropDown            matlab.ui.control.DropDown
        CombinationDropDownLabel       matlab.ui.control.Label
        RebuildImagesButton            matlab.ui.control.Button
        WindowSizeEditField            matlab.ui.control.NumericEditField
        WindowSizeEditFieldLabel       matlab.ui.control.Label
        WindowDropDown                 matlab.ui.control.DropDown
        WindowDropDownLabel            matlab.ui.control.Label
        FFTDropDown                    matlab.ui.control.DropDown
        Sizeof2DFFTLabel               matlab.ui.control.Label
        axis_coordinates               matlab.ui.control.UIAxes
        axis_image                     matlab.ui.control.UIAxes
        axis_kspace                    matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        currentscan % Currently selected scan
        scan_data %actual processed scan data
        currentmap
        mask

        ImagePlaneKey = []
        ImageHandle = []
        KspaceHandle = []
        OverlayHandle = []
        CompassKey = []
        KspaceKey = []
        ViewStates = {}
        CurrentTaskKind = ''
        JobRunner = []
        FilterDialogs = []
        DatasetGeneration = uint64(0)
        BusyControls = {}
        svsOverlayScan = [] % most recently selected SVS prescription
        svsLocalizerScan = [] % imaging scan selected from Scan List

        % cvi42 ventricular-contour overlay.  Source contours remain in
        % clinical patient coordinates; registration is applied only when
        % drawing them on an FCI image.
        cviData = []
        cviRegistration = struct()
        CVIOverlayHandles = []
        cviDialog = []
        cviMenuHandles = struct()
    end
    
    methods (Access = private)
        
        function   activeimage(app)
            app.currentscan = app.ScanList.Value;
            % The viewer consumes combined images; the public core API can
            % still retain per-receiver RSSQ images for external analyses.
            app.scan_data{app.currentscan}.retain_coil_images = false;
            if scanIsSVS(app, app.scan_data{app.currentscan})
                app.svsOverlayScan = app.currentscan;
            elseif scanIsImaging(app, app.scan_data{app.currentscan})
                % Selecting an imaging entry explicitly makes that dataset
                % the localizer. The most recently selected SVS prescription
                % remains available as an overlay on every imaging dataset.
                app.svsLocalizerScan = app.currentscan;
            end
            if isempty(app.scan_data{app.currentscan}.rot_theta)
                app.scan_data{app.currentscan}.rot_theta = 0;
            end

            if ~isprop(app.scan_data{app.currentscan}, 'geomT') || ...
                    isempty(app.scan_data{app.currentscan}.geomT)
              app.scan_data{app.currentscan}.geomT = eye(3);
            end

            %generate string lists with appropriate evolution fields, times etc
            fft_list = cellfun(@num2str,{app.scan_data{app.currentscan}.samples,app.scan_data{app.currentscan}.samples*2,app.scan_data{app.currentscan}.samples*4},'UniformOutput',false);
            bevol_list = sprintfc('%g',round(app.scan_data{app.currentscan}.fieldpoints,2));
            tevol_list = sprintfc('%g',round(app.scan_data{app.currentscan}.timepoints(1,:),2));
            
            
            %populate the GUI with the lists
            app.FFTDropDown.Items = fft_list;
            app.EvolutionFieldDropDown.Items = bevol_list;
            app.EvolutionTimeDropDown.Items = tevol_list;
            app.SliceDropDown.Items = sprintfc('%i',1:app.scan_data{app.currentscan}.slices);
            
            %allows use of indexing
            app.EvolutionFieldDropDown.ItemsData = 1:numel(app.EvolutionFieldDropDown.Items); %allows us to get indexes
            app.EvolutionTimeDropDown.ItemsData = 1:numel(app.EvolutionTimeDropDown.Items);
            app.SliceDropDown.ItemsData = 1:app.scan_data{app.currentscan}.slices;
            % Scanner receiver IDs are zero-based, whereas MATLAB array
            % indices are one-based.  Display 0..N-1 but retain 1..N as
            % ItemsData so callbacks can index the receiver dimension safely.
            app.MultichannelDropDown.Items = string(0:(app.scan_data{app.currentscan}.n_receivers-1));
            app.MultichannelDropDown.ItemsData = 1:numel(app.MultichannelDropDown.Items);
            
            if length(app.SliceDropDown.ItemsData) ==1
                app.SliceDropDown.Enable = 0;
            else
                app.SliceDropDown.Enable = 1;
            end
            
            if length(app.EvolutionTimeDropDown.ItemsData) ==1
                app.EvolutionTimeDropDown.Enable = 0;
                app.GenerateROIDispersionCurveButton.Enable=0;
                app.GenerateRelaxationMapsButton.Enable =0;
                
                
            else
                app.EvolutionTimeDropDown.Enable = 1;
                app.GenerateROIDispersionCurveButton.Enable=1;
                app.GenerateRelaxationMapsButton.Enable =1;
                
            end
            
            if length(app.EvolutionFieldDropDown.ItemsData) ==1
                app.EvolutionFieldDropDown.Enable = 0;
            else
                app.EvolutionFieldDropDown.Enable = 1;
            end
            
            app.FFTDropDown.Value =  num2str(app.scan_data{app.currentscan}.fft_size);
            update3DFFTDropDown(app);
            app.WindowSizeEditField.Value = app.scan_data{app.currentscan}.window_size;
            
            if strcmp(app.scan_data{app.currentscan}.window_function, 'none')
                app.scan_data{app.currentscan}.window_function = 'None';
            end
            app.WindowDropDown.Value = app.scan_data{app.currentscan}.window_function;
            
            if isempty(app.scan_data{app.currentscan}.denoise_filter)
                app.scan_data{app.currentscan}.denoise_filter = 'None';
            end

            if isempty(app.scan_data{app.currentscan}.denoise_params) || ...
                    ~isstruct(app.scan_data{app.currentscan}.denoise_params)
                app.scan_data{app.currentscan}.denoise_params = ...
                    ffc_default_denoise_params(app.scan_data{app.currentscan}.denoise_filter);
            end

            app.FilterSettingsButton.Enable = 'on';
            app.FilterSettingsButton.Tooltip = {sprintf('Current denoise filter: %s', ...
                app.scan_data{app.currentscan}.denoise_filter)};

            nCh = app.scan_data{app.currentscan}.n_receivers;
            selectedReceivers = app.scan_data{app.currentscan}.multichannel_recon;
            if isempty(selectedReceivers)
                selectedReceivers = 1:nCh;
            else
                if islogical(selectedReceivers), selectedReceivers = find(selectedReceivers(:).'); end
                selectedReceivers = unique(round(double(selectedReceivers(:).')), 'stable');
                selectedReceivers = selectedReceivers(isfinite(selectedReceivers) & ...
                    selectedReceivers >= 1 & selectedReceivers <= nCh);
                if isempty(selectedReceivers)
                    selectedReceivers = 1:nCh;
                end
            end
            app.scan_data{app.currentscan}.multichannel_recon = selectedReceivers;

            if isempty(app.scan_data{app.currentscan}.biascorrect) || ...
                    ~ismember(string(app.scan_data{app.currentscan}.biascorrect), ["Off", "ACS"])
                app.scan_data{app.currentscan}.biascorrect = 'Off';
            end
            app.BiasFieldDropDown.Value = app.scan_data{app.currentscan}.biascorrect;

            if isempty(app.scan_data{app.currentscan}.combination) || ...
                    ~ismember(string(app.scan_data{app.currentscan}.combination), ...
                    ["RSSQ", "Adaptive", "ACS/Walsh"])
                app.scan_data{app.currentscan}.combination = 'Adaptive';
            end
            app.CombinationDropDown.Value = app.scan_data{app.currentscan}.combination;

            whiteningCheckBox = getNoiseWhiteningCheckBox(app);
            if nCh <= 1
                app.scan_data{app.currentscan}.noise_whitening = false;
                whiteningCheckBox.Value = false;
                whiteningCheckBox.Enable = 'off';
                whiteningCheckBox.Tooltip = {'Noise whitening requires at least two acquired receive channels.'};
                app.CombinationDropDown.Enable = 'off';
                app.CombinationDropDownLabel.Enable = 'off';
            else
                if isempty(app.scan_data{app.currentscan}.noise_whitening)
                    app.scan_data{app.currentscan}.noise_whitening = true;
                end
                whiteningCheckBox.Value = logical(app.scan_data{app.currentscan}.noise_whitening);
                whiteningCheckBox.Enable = 'on';
                whiteningCheckBox.Tooltip = {'Prewhiten all acquired receivers, then reconstruct only the selected receiver IDs.'};
                app.CombinationDropDown.Enable = 'on';
                app.CombinationDropDownLabel.Enable = 'on';
            end

            app.SliceDropDown.Value = 1;
            app.EvolutionFieldDropDown.Value =1;
            app.EvolutionTimeDropDown.Value =1;
            app.MultichannelDropDown.Items     = string(0:(nCh-1));
            app.MultichannelDropDown.ItemsData = 1:nCh;
            app.MultichannelDropDown.Value     = selectedReceivers;
            %kspace to image space with selected recon parameters
            if ~ffc_reconstruction_valid(app.scan_data{app.currentscan})
                startScanTask(app,'reconstruct',struct('force',false));
                return;
            end
            configureViewOrientation(app);
            updateSliceDropDownFromImages(app);
            [img0, ~, hasLocalizer] = currentMagnitudePlane(app);
            if ~restoreViewState(app) && ~isempty(img0)
                autoWindowFromImage(app,img0);
            end

            if scanIsSVS(app, app.scan_data{app.currentscan})
                app.GenerateRelaxationMapsButton.Enable = 'off';
                app.ImageTypeDropDown.Value = 'Magnitude';
                app.ImageTypeDropDown.Enable = 'off';
                app.FFTDropDown.Enable = 'off';
                app.FFTDropDown_3D.Enable = 'off';
                app.WindowDropDown.Enable = 'off';
                app.WindowSizeEditField.Enable = 'off';
                app.RegisterImagesButton.Enable = 'off';
                if ~hasLocalizer
                    app.SliceDropDown.Enable = 'off';
                end
            else
                app.ImageTypeDropDown.Enable = 'on';
                app.FFTDropDown.Enable = 'on';
                app.WindowDropDown.Enable = 'on';
                app.WindowSizeEditField.Enable = 'on';
                app.RegisterImagesButton.Enable = 'on';
                if app.scan_data{app.currentscan}.n_timepoints > 1
                    app.GenerateRelaxationMapsButton.Enable = 'on';
                end
                update3DFFTDropDown(app);
            end
            
            drawimages(app)



        end

        function restoreMouseHandlers(app,f,down,move,up,wheel)
            if ~isgraphics(f), return; end
            f.WindowButtonDownFcn = down; f.WindowButtonMotionFcn = move;
            f.WindowButtonUpFcn = up; f.WindowScrollWheelFcn = wheel;
        end

        function limits = windowLimits(app)
            level = app.WindowLevelSlider.Value;
            width = max(app.WindowWidthSlider.Value,eps(max(1,abs(level))));
            limits = level + [-0.5 0.5]*width;
        end

        function updateWindowOnly(app)
            if isempty(app.ImageHandle) || (isempty(app.ImageHandle) || ~isgraphics(app.ImageHandle)), return; end
            if strcmp(app.ImageTypeDropDown.Value,'Magnitude')
                app.axis_image.CLim = windowLimits(app);
                drawnow limitrate;
            end
        end

        function mainWindowChanging(app,source,event)
            source.Value = event.Value;
            updateWindowOnly(app);
        end

        function saveViewState(app)
            if isempty(app.currentscan) || isempty(app.scan_data), return; end
            app.ViewStates{app.currentscan} = struct('slice',app.SliceDropDown.Value, ...
                'time',app.EvolutionTimeDropDown.Value,'field',app.EvolutionFieldDropDown.Value, ...
                'level',app.WindowLevelSlider.Value,'width',app.WindowWidthSlider.Value, ...
                'levelLimits',app.WindowLevelSlider.Limits,'widthLimits',app.WindowWidthSlider.Limits);
        end

        function restored = restoreViewState(app)
            restored = false;
            if numel(app.ViewStates) < app.currentscan || isempty(app.ViewStates{app.currentscan}), return; end
            state = app.ViewStates{app.currentscan};
            app.EvolutionFieldDropDown.Value = min(state.field,max(app.EvolutionFieldDropDown.ItemsData));
            scan = app.scan_data{app.currentscan};
            times = scan.timepoints(min(app.EvolutionFieldDropDown.Value,size(scan.timepoints,1)),:);
            app.EvolutionTimeDropDown.Items = sprintfc('%g',round(times,2));
            app.EvolutionTimeDropDown.ItemsData = 1:numel(times);
            app.EvolutionTimeDropDown.Value = min(state.time,numel(times));
            app.SliceDropDown.Value = min(state.slice,max(app.SliceDropDown.ItemsData));
            app.WindowLevelSlider.Limits = state.levelLimits;
            app.WindowWidthSlider.Limits = state.widthLimits;
            app.WindowLevelSlider.Value = state.level;
            app.WindowWidthSlider.Value = state.width;
            restored = true;
        end

        function setTaskBusy(app,busy,label)
            if nargin < 3, label = ''; end
            if busy
                controls = {app.FFTDropDown,app.FFTDropDown_3D,app.WindowDropDown, ...
                    app.WindowSizeEditField,app.MultichannelDropDown,app.CombinationDropDown, ...
                    app.BiasFieldDropDown,getNoiseWhiteningCheckBox(app), ...
                    app.FilterSettingsButton,app.RegisterImagesButton, ...
                    app.GenerateRelaxationMapsButton,app.GenerateROIDispersionCurveButton};
                app.BusyControls = cell(numel(controls),2);
                for i = 1:numel(controls)
                    app.BusyControls(i,:) = {controls{i},controls{i}.Enable};
                    controls{i}.Enable = 'off';
                end
                app.RebuildImagesButton.Text = 'Cancel';
                app.UIFigure.Name = ['FFC Processing - ' label '...'];
            else
                for i = 1:size(app.BusyControls,1)
                    if isvalid(app.BusyControls{i,1})
                        app.BusyControls{i,1}.Enable = app.BusyControls{i,2};
                    end
                end
                app.BusyControls = {};
                app.RebuildImagesButton.Text = 'Rebuild Images';
                app.UIFigure.Name = 'FFC Processing - Performance Edition + CVI contours TEST 2';
            end
        end

        function cancelScanTask(app)
            if ~isempty(app.JobRunner) && app.JobRunner.isRunning()
                rollbackReconstructionSettings(app);
                app.JobRunner.cancel();
            end
            app.CurrentTaskKind = '';
            setTaskBusy(app,false);
        end

        function rollbackReconstructionSettings(app)
            if ~strcmp(app.CurrentTaskKind,'reconstruct') || isempty(app.currentscan), return; end
            scan = app.scan_data{app.currentscan};
            if isempty(scan.reconstruction_key), return; end
            names = {'fft_size','fft_size_3d','window_size','window_function', ...
                'multichannel_recon','combination','biascorrect','noise_whitening'};
            for i = 1:numel(names), scan.(names{i}) = scan.reconstruction_key.(names{i}); end
            app.scan_data{app.currentscan} = scan;
        end

        function closeFilterDialogs(app)
            for i = 1:numel(app.FilterDialogs)
                if isgraphics(app.FilterDialogs(i)), close(app.FilterDialogs(i)); end
            end
            app.FilterDialogs = [];
        end

        function startScanTask(app,kind,payload)
            cancelScanTask(app);
            if ~strcmp(kind,'filter'), closeFilterDialogs(app); end
            app.CurrentTaskKind = kind;
            index = app.currentscan;
            dataset = app.DatasetGeneration;
            sourceGeneration = app.scan_data{index}.magnitude_generation;
            if strcmp(kind,'reconstruct') && isempty(app.scan_data{index}.magimage)
                cla(app.axis_image); cla(app.axis_kspace);
                text(app.axis_image,0.5,0.5,'Reconstructing...', ...
                    'Units','normalized','HorizontalAlignment','center');
            end
            setTaskBusy(app,true,kind);
            app.JobRunner.start(kind,app.scan_data{index},payload, ...
                @(result) finishScanTask(app,result,index,dataset,sourceGeneration,kind), ...
                @(err) failScanTask(app,err));
        end

        function finishScanTask(app,result,index,dataset,sourceGeneration,kind)
            if ~isvalid(app) || ~isgraphics(app.UIFigure), return; end
            app.CurrentTaskKind = '';
            setTaskBusy(app,false);
            if dataset ~= app.DatasetGeneration || index > numel(app.scan_data) || ...
                    sourceGeneration ~= app.scan_data{index}.magnitude_generation
                return;
            end
            % Display transforms may change while numerical work is running.
            for name = {'geomT','rot_theta','flipH','flipV','view_orientation','R1T1'}
                result.(name{1}) = app.scan_data{index}.(name{1});
            end
            app.scan_data{index} = result;
            if index == app.currentscan
                if strcmp(kind,'reconstruct')
                    activeimage(app);
                else
                    drawimages(app);
                    updateFilterButtonTooltip(app);
                end
            end
        end

        function failScanTask(app,err)
            if ~isvalid(app) || ~isgraphics(app.UIFigure), return; end
            rollbackReconstructionSettings(app);
            app.CurrentTaskKind = '';
            setTaskBusy(app,false);
            if ~isempty(app.currentscan) && ffc_reconstruction_valid(app.scan_data{app.currentscan})
                activeimage(app);
            end
            fprintf(2,'%s\n',getReport(err,'extended','hyperlinks','off'));
            uialert(app.UIFigure,err.message,'Processing failed');
        end


        function checkBox = getNoiseWhiteningCheckBox(app)
            if ~isappdata(app.panel_recon, 'NoiseWhiteningCheckBox')
                error('FFCProcessing:MissingWhiteningControl', ...
                    'The noise-whitening control has not been initialised.');
            end
            checkBox = getappdata(app.panel_recon, 'NoiseWhiteningCheckBox');
        end

        function dropDown = getViewOrientationDropDown(app)
            if ~isappdata(app.CurrentImagePanel, 'ViewOrientationDropDown')
                error('FFCProcessing:MissingOrientationControl', ...
                    'The anatomical-view control has not been initialised.');
            end
            dropDown = getappdata(app.CurrentImagePanel, ...
                'ViewOrientationDropDown');
        end

        function tf = scanIsSVS(app, scan) %#ok<INUSL>
            tf = false;
            try
                if isprop(scan, 'is_svs') && ~isempty(scan.is_svs)
                    tf = logical(scan.is_svs);
                end
            catch
            end
            if tf
                return;
            end
            try
                tf = contains(lower(string(scan.sequence)), 'svs');
            catch
            end
        end

        function tf = scanIsImaging(app, scan)
            tf = false;
            if scanIsSVS(app, scan)
                return;
            end
            try
                tf = scan.samples > 1 && scan.views > 1;
            catch
            end
        end

        function index = findSVSOverlayScan(app)
            index = 0;
            if isempty(app.scan_data)
                return;
            end
            if ~isempty(app.currentscan) && ...
                    scanIsSVS(app, app.scan_data{app.currentscan})
                index = app.currentscan;
                return;
            end
            candidate = app.svsOverlayScan;
            if ~isempty(candidate) && isscalar(candidate) && ...
                    candidate >= 1 && candidate <= numel(app.scan_data) && ...
                    scanIsSVS(app, app.scan_data{candidate})
                index = candidate;
            end
        end

        function index = findSVSLocalizer(app)
            index = 0;
            if isempty(app.scan_data)
                return;
            end
            selected = app.svsLocalizerScan;
            if ~isempty(selected) && isscalar(selected) && ...
                    selected >= 1 && selected <= numel(app.scan_data) && ...
                    scanIsImaging(app, app.scan_data{selected})
                index = selected;
                return;
            end
            candidates = [];
            true3D = [];
            for scanIndex = 1:numel(app.scan_data)
                scan = app.scan_data{scanIndex};
                if ~scanIsImaging(app, scan)
                    continue;
                end
                candidates(end + 1) = scanIndex; %#ok<AGROW>
                try
                    if scan.TwoDimensional ~= 1 && scan.slices > 1
                        true3D(end + 1) = scanIndex; %#ok<AGROW>
                    end
                catch
                end
            end
            if ~isempty(true3D)
                index = true3D(1);
            elseif ~isempty(candidates)
                index = candidates(1);
            end
            app.svsLocalizerScan = index;
        end

        function prepareLocalizer(app, index)
            if index < 1 || index > numel(app.scan_data)
                return;
            end
            localizer = app.scan_data{index};
            if ~isempty(localizer.magimage)
                return;
            end
            nChannels = max(1, double(localizer.n_receivers));
            if isempty(localizer.multichannel_recon)
                localizer.multichannel_recon = 1:nChannels;
            end
            if isempty(localizer.combination)
                localizer.combination = 'Adaptive';
            end
            if isempty(localizer.biascorrect)
                localizer.biascorrect = 'Off';
            end
            if isempty(localizer.noise_whitening)
                localizer.noise_whitening = nChannels > 1;
            end
            localizer.retain_coil_images = false;
            localizer = buildimages(localizer,false);
            app.scan_data{index} = localizer;
        end

        function [displayScan, displayIndex, isLocalizer] = getDisplayScan(app)
            displayIndex = app.currentscan;
            isLocalizer = false;
            currentScan = app.scan_data{app.currentscan};
            if scanIsSVS(app, currentScan)
                localizerIndex = findSVSLocalizer(app);
                if localizerIndex > 0
                    prepareLocalizer(app, localizerIndex);
                    displayIndex = localizerIndex;
                    isLocalizer = true;
                end
            end
            displayScan = app.scan_data{displayIndex};
        end

        function planeName = currentViewOrientation(app)
            dropDown = getViewOrientationDropDown(app);
            planeName = char(string(dropDown.Value));
        end

        function geomT = currentDisplayTransform(app)
            [displayScan, ~, ~] = getDisplayScan(app);
            geomT = displayScan.geomT;
            if isempty(geomT)
                geomT = eye(3);
            end
        end

        function planeName = acquiredPlaneName(app, scan) %#ok<INUSL>
            planeName = 'Axial';
            try
                orientation = double(scan.orientation(:));
                if numel(orientation) ~= 6
                    return;
                end
                normal = cross(orientation(1:3), orientation(4:6));
                [~, dimension] = max(abs(normal));
                names = {'Sagittal', 'Coronal', 'Axial'};
                planeName = names{dimension};
            catch
            end
        end

        function configureViewOrientation(app)
            dropDown = getViewOrientationDropDown(app);
            label = getappdata(app.CurrentImagePanel, 'ViewOrientationLabel');
            [displayScan, ~, ~] = getDisplayScan(app);
            geometry = ffc_scan_spatial_geometry(displayScan);
            currentScan = app.scan_data{app.currentscan};
            value = '';
            try
                value = char(string(currentScan.view_orientation));
            catch
            end
            if ~ismember(value, {'Axial', 'Coronal', 'Sagittal'})
                value = acquiredPlaneName(app, displayScan);
            end
            if ~geometry.is3D
                value = acquiredPlaneName(app, displayScan);
                dropDown.Enable = 'off';
                dropDown.Tooltip = ...
                    {'Orthogonal reformatting requires a true 3D acquisition.'};
                label.Enable = 'off';
            else
                dropDown.Enable = 'on';
                dropDown.Tooltip = ...
                    {'Display the reconstructed 3D volume in an anatomical plane.'};
                label.Enable = 'on';
            end
            dropDown.Value = value;
            currentScan.view_orientation = value;
            app.scan_data{app.currentscan} = currentScan;
        end

        function [image, view, hasLocalizer] = currentMagnitudePlane(app)
            [displayScan, displayIndex, hasLocalizer] = getDisplayScan(app);
            if scanIsSVS(app, app.scan_data{app.currentscan}) && ...
                    displayIndex == app.currentscan
                image = [];
                view = ffc_anatomical_plane_geometry(displayScan, ...
                    acquiredPlaneName(app, displayScan), 1, ...
                    currentDisplayTransform(app));
                return;
            end
            if hasLocalizer
                timeIndex = 1;
                fieldIndex = 1;
            else
                timeIndex = min(app.EvolutionTimeDropDown.Value, ...
                    size(displayScan.magimage, 4));
                fieldIndex = min(app.EvolutionFieldDropDown.Value, ...
                    size(displayScan.magimage, 5));
            end
            volume = displayScan.magimage(:, :, :, timeIndex, fieldIndex);
            [image, view] = ffc_extract_anatomical_plane(volume, displayScan, ...
                currentViewOrientation(app), app.SliceDropDown.Value, ...
                currentDisplayTransform(app));
        end

        function drawOrientationCompass(app, view)
            compassKey = [view.basisXPatient,view.basisYPatient,view.basisNormalPatient];
            if isequaln(app.CompassKey,compassKey) && ~isempty(app.axis_coordinates.Children), return; end
            app.CompassKey = compassKey;
            cla(app.axis_coordinates);
            hold(app.axis_coordinates, 'on');
            axis(app.axis_coordinates, [-1.2 1.2 -1.2 1.2 -1.2 1.2]);
            axis(app.axis_coordinates, 'off');
            grid(app.axis_coordinates, 'off');

            % CData rows increase downward in the image axes.  Use the
            % corresponding screen-up direction for the compass so that the
            % labels describe the pixels the viewer sees, not the array index.
            patientToDisplay = [view.basisXPatient, ...
                -view.basisYPatient, view.basisNormalPatient].';
            % RS2D subject convention: +X=R, +Y=A and +Z=H (therefore
            % -X=L, -Y=P and -Z=F).
            L = patientToDisplay * [-1; 0; 0];
            R = -L;
            A = patientToDisplay * [0; 1; 0];
            P = -A;
            H = patientToDisplay * [0; 0; 1];
            F = -H;

            plot3(app.axis_coordinates, [L(1) R(1)], [L(2) R(2)], ...
                [L(3) R(3)], 'LineWidth', 3);
            plot3(app.axis_coordinates, [P(1) A(1)], [P(2) A(2)], ...
                [P(3) A(3)], 'LineWidth', 3);
            plot3(app.axis_coordinates, [F(1) H(1)], [F(2) H(2)], ...
                [F(3) H(3)], 'LineWidth', 3);

            scale = 1.15;
            labels = {'L', 'R', 'P', 'A', 'H', 'F'};
            points = [L, R, P, A, H, F];
            for labelIndex = 1:numel(labels)
                text(app.axis_coordinates, ...
                    points(1, labelIndex) .* scale, ...
                    points(2, labelIndex) .* scale, ...
                    points(3, labelIndex) .* scale, labels{labelIndex}, ...
                    'FontWeight', 'bold', 'FontSize', 14);
            end
            hold(app.axis_coordinates, 'off');
        end

        function drawSVSVoxelOverlay(app, view)
            if (isempty(app.OverlayHandle) || ~isgraphics(app.OverlayHandle))
                app.OverlayHandle = line(app.axis_image,nan,nan,'Color','r', ...
                    'LineWidth',2.5,'HitTest','off','PickableParts','none');
            end
            app.OverlayHandle.Visible = 'off';
            svsIndex = findSVSOverlayScan(app);
            if svsIndex == 0, return; end
            [x,y,details] = ffc_voxel_plane_intersection(app.scan_data{svsIndex},view);
            if details.intersects
                set(app.OverlayHandle,'XData',x,'YData',y,'Visible','on');
            end
        end

        function drawCVIContourOverlay(app, view)
            if isempty(app.CVIOverlayHandles) || ...
                    numel(app.CVIOverlayHandles) ~= 3 || ...
                    any(~isgraphics(app.CVIOverlayHandles))
                app.CVIOverlayHandles = gobjects(1,3);
                app.CVIOverlayHandles(1) = line(app.axis_image,nan,nan, ...
                    'Color',[1 0 0.45],'LineWidth',2.0, ...
                    'HitTest','off','PickableParts','none');
                app.CVIOverlayHandles(2) = line(app.axis_image,nan,nan, ...
                    'Color',[0 0.9 0.9],'LineWidth',2.0, ...
                    'HitTest','off','PickableParts','none');
                app.CVIOverlayHandles(3) = line(app.axis_image,nan,nan, ...
                    'Color',[0.95 0.85 0],'LineWidth',2.0, ...
                    'HitTest','off','PickableParts','none');
            end
            set(app.CVIOverlayHandles,'Visible','off');
            if isempty(app.cviData) || isempty(fieldnames(app.cviRegistration)) || ...
                    ~isfield(app.cviRegistration,'visible') || ...
                    ~app.cviRegistration.visible
                return;
            end
            seriesIndex = min(max(round(app.cviRegistration.seriesIndex),1), ...
                numel(app.cviData.series));
            series = app.cviData.series(seriesIndex);
            phaseIndex = app.cviRegistration.phaseIndex;
            types = {'LVEndo','LVEpi','RVEndo'};
            show = [app.cviRegistration.showLVEndo, ...
                app.cviRegistration.showLVEpi,app.cviRegistration.showRVEndo];
            for t = 1:numel(types)
                if ~show(t), continue; end
                matches = strcmp({series.Surfaces.Type},types{t}) & ...
                    [series.Surfaces.PhaseIndex] == phaseIndex;
                surfaceIndex = find(matches,1);
                if isempty(surfaceIndex), continue; end
                [x,y] = cvi_surface_plane_intersection( ...
                    series.Surfaces(surfaceIndex),view,app.cviRegistration, ...
                    series.SourceCenterPatientMm);
                if ~isempty(x)
                    set(app.CVIOverlayHandles(t),'XData',x,'YData',y,'Visible','on');
                end
            end
        end

        function importCVIContours(app)
            [file,path] = uigetfile( ...
                {'*.cvi42ws;*.cvi42wsx','cvi42 workspace (*.cvi42ws, *.cvi42wsx)'}, ...
                'Select cvi42 workspace');
            if isequal(file,0), return; end
            dicomRoot = uigetdir(path, ...
                'Select clinical DICOM folder used by this cvi42 workspace');
            if isequal(dicomRoot,0), return; end
            progress = uiprogressdlg(app.UIFigure,'Title','Importing cvi42 contours', ...
                'Message','Decoding workspace and matching DICOM geometry...', ...
                'Indeterminate','on');
            cleaner = onCleanup(@() cvi_close_progress(progress)); %#ok<NASGU>
            try
                cvi = cvi_import_contours(fullfile(path,file),dicomRoot);
            catch err
                fprintf(2,'%s\n',getReport(err,'extended','hyperlinks','off'));
                uialert(app.UIFigure,err.message,'CVI contour import failed');
                return;
            end
            app.cviData = cvi;
            seriesIndex = cvi.defaultSeriesIndex;
            series = cvi.series(seriesIndex);
            phaseIndex = series.PhaseIndices(1);
            ed = find(startsWith(string(series.PhaseLabels),'ED'),1);
            if ~isempty(ed), phaseIndex = series.PhaseIndices(ed); end
            app.cviRegistration = struct( ...
                'visible',true,'seriesIndex',seriesIndex,'phaseIndex',phaseIndex, ...
                'translationMm',[0 0 0],'rotationDeg',[0 0 0], ...
                'showLVEndo',true,'showLVEpi',true,'showRVEndo',true);
            if isfield(app.cviMenuHandles,'show') && isgraphics(app.cviMenuHandles.show)
                app.cviMenuHandles.show.Checked = 'on';
            end
            if ~isempty(app.scan_data) && ~isempty(app.currentscan)
                drawimages(app);
            end
            openCVIRegistrationDialog(app);
        end

        function toggleCVIContours(app)
            if isempty(app.cviData)
                importCVIContours(app);
                return;
            end
            app.cviRegistration.visible = ~app.cviRegistration.visible;
            if isfield(app.cviMenuHandles,'show') && isgraphics(app.cviMenuHandles.show)
                if app.cviRegistration.visible
                    app.cviMenuHandles.show.Checked = 'on';
                else
                    app.cviMenuHandles.show.Checked = 'off';
                end
            end
            if ~isempty(app.scan_data) && ~isempty(app.currentscan), drawimages(app); end
        end

        function clearCVIContours(app)
            app.cviData = [];
            app.cviRegistration = struct();
            if ~isempty(app.CVIOverlayHandles)
                delete(app.CVIOverlayHandles(isgraphics(app.CVIOverlayHandles)));
            end
            app.CVIOverlayHandles = [];
            if ~isempty(app.cviDialog) && isgraphics(app.cviDialog), delete(app.cviDialog); end
            app.cviDialog = [];
            if isfield(app.cviMenuHandles,'show') && isgraphics(app.cviMenuHandles.show)
                app.cviMenuHandles.show.Checked = 'off';
            end
            if ~isempty(app.scan_data) && ~isempty(app.currentscan), drawimages(app); end
        end

        function openCVIRegistrationDialog(app)
            if isempty(app.cviData)
                uialert(app.UIFigure,'Import a cvi42 workspace first.','No CVI contours');
                return;
            end
            if ~isempty(app.cviDialog) && isgraphics(app.cviDialog)
                app.cviDialog.Visible = 'on';
                updateCVIDialog(app);
                return;
            end
            d = uifigure('Name','CVI Contour Registration', ...
                'Position',[180 180 430 455],'Resize','off');
            app.cviDialog = d;

            uilabel(d,'Text','Clinical series','Position',[20 410 95 22]);
            seriesDrop = uidropdown(d,'Position',[120 410 285 22], ...
                'ValueChangedFcn',@(src,event) cviSeriesChanged(app,src));
            setappdata(d,'SeriesDrop',seriesDrop);
            uilabel(d,'Text','Cardiac phase','Position',[20 375 95 22]);
            phaseDrop = uidropdown(d,'Position',[120 375 285 22], ...
                'ValueChangedFcn',@(src,event) cviPhaseChanged(app,src));
            setappdata(d,'PhaseDrop',phaseDrop);

            showEndo = uicheckbox(d,'Text','LV endocardium','Position',[20 335 120 22], ...
                'ValueChangedFcn',@(src,event) cviToggleChanged(app,'showLVEndo',src.Value));
            showEpi = uicheckbox(d,'Text','LV epicardium','Position',[150 335 120 22], ...
                'ValueChangedFcn',@(src,event) cviToggleChanged(app,'showLVEpi',src.Value));
            showRV = uicheckbox(d,'Text','RV endocardium','Position',[280 335 125 22], ...
                'ValueChangedFcn',@(src,event) cviToggleChanged(app,'showRVEndo',src.Value));
            setappdata(d,'ShowEndo',showEndo); setappdata(d,'ShowEpi',showEpi); setappdata(d,'ShowRV',showRV);

            uilabel(d,'Text','Rigid registration (clinical contour stack -> FCI)', ...
                'FontWeight','bold','Position',[20 298 360 22]);
            labels = {'X translation (mm)','Y translation (mm)','Z translation (mm)', ...
                'X rotation (deg)','Y rotation (deg)','Z rotation (deg)'};
            fields = gobjects(1,6);
            for k = 1:6
                y = 298-k*34;
                uilabel(d,'Text',labels{k},'Position',[35 y 145 22]);
                fields(k) = uieditfield(d,'numeric','Position',[190 y 95 22], ...
                    'ValueChangedFcn',@(src,event) cviNumericChanged(app,k,src.Value));
            end
            setappdata(d,'NumericFields',fields);

            uibutton(d,'Text','Reset registration','Position',[300 145 105 28], ...
                'ButtonPushedFcn',@(src,event) resetCVIRegistration(app));
            visible = uicheckbox(d,'Text','Show contours','Position',[20 110 120 22], ...
                'ValueChangedFcn',@(src,event) cviVisibilityChanged(app,src.Value));
            setappdata(d,'VisibleCheck',visible);
            status = uilabel(d,'Text','', 'Position',[20 48 385 52], ...
                'WordWrap','on');
            setappdata(d,'StatusLabel',status);
            uibutton(d,'Text','Close','Position',[315 15 90 28], ...
                'ButtonPushedFcn',@(src,event) delete(d));
            updateCVIDialog(app);
        end

        function updateCVIDialog(app)
            if isempty(app.cviDialog) || ~isgraphics(app.cviDialog) || isempty(app.cviData), return; end
            d = app.cviDialog;
            seriesDrop = getappdata(d,'SeriesDrop');
            items = arrayfun(@(s) sprintf('%s (%d contours)',s.Description,s.ContourCount), ...
                app.cviData.series,'UniformOutput',false);
            seriesDrop.Items = items;
            seriesDrop.ItemsData = 1:numel(items);
            seriesDrop.Value = app.cviRegistration.seriesIndex;
            series = app.cviData.series(app.cviRegistration.seriesIndex);
            phaseDrop = getappdata(d,'PhaseDrop');
            phaseDrop.Items = series.PhaseLabels;
            phaseDrop.ItemsData = series.PhaseItemsData;
            if any(series.PhaseItemsData == app.cviRegistration.phaseIndex)
                phaseDrop.Value = app.cviRegistration.phaseIndex;
            else
                app.cviRegistration.phaseIndex = series.PhaseItemsData(1);
                phaseDrop.Value = app.cviRegistration.phaseIndex;
            end
            getappdata(d,'ShowEndo').Value = app.cviRegistration.showLVEndo;
            getappdata(d,'ShowEpi').Value = app.cviRegistration.showLVEpi;
            getappdata(d,'ShowRV').Value = app.cviRegistration.showRVEndo;
            fields = getappdata(d,'NumericFields');
            values = [app.cviRegistration.translationMm,app.cviRegistration.rotationDeg];
            for k=1:6, fields(k).Value = values(k); end
            getappdata(d,'VisibleCheck').Value = app.cviRegistration.visible;
            status = getappdata(d,'StatusLabel');
            status.Text = sprintf(['Imported %d contours. Selected series: %d frames, %d cardiac phases. ', ...
                '%d workspace contour(s) could not be matched to DICOM geometry.'], ...
                numel(app.cviData.contours),series.NumberOfFrames,series.NumberOfPhases, ...
                app.cviData.unmatchedContourCount);
        end

        function cviSeriesChanged(app,source)
            app.cviRegistration.seriesIndex = source.Value;
            series = app.cviData.series(source.Value);
            app.cviRegistration.phaseIndex = series.PhaseItemsData(1);
            ed = find(startsWith(string(series.PhaseLabels),'ED'),1);
            if ~isempty(ed), app.cviRegistration.phaseIndex = series.PhaseItemsData(ed); end
            updateCVIDialog(app); redrawCVI(app);
        end

        function cviPhaseChanged(app,source)
            app.cviRegistration.phaseIndex = source.Value;
            redrawCVI(app);
        end

        function cviToggleChanged(app,field,value)
            app.cviRegistration.(field) = logical(value);
            redrawCVI(app);
        end

        function cviVisibilityChanged(app,value)
            app.cviRegistration.visible = logical(value);
            if isfield(app.cviMenuHandles,'show') && isgraphics(app.cviMenuHandles.show)
                if value, app.cviMenuHandles.show.Checked='on'; else, app.cviMenuHandles.show.Checked='off'; end
            end
            redrawCVI(app);
        end

        function cviNumericChanged(app,index,value)
            if index <= 3
                app.cviRegistration.translationMm(index) = value;
            else
                app.cviRegistration.rotationDeg(index-3) = value;
            end
            redrawCVI(app);
        end

        function resetCVIRegistration(app)
            app.cviRegistration.translationMm = [0 0 0];
            app.cviRegistration.rotationDeg = [0 0 0];
            updateCVIDialog(app); redrawCVI(app);
        end

        function redrawCVI(app)
            if ~isempty(app.scan_data) && ~isempty(app.currentscan), drawimages(app); end
        end

        function applyImageFovGeometry(app, view)
            imageHandles = app.ImageHandle;
            if isempty(imageHandles) || isempty(view) || ...
                    ~isfield(view, 'fovX') || ~isfield(view, 'fovY')
                return;
            end
            fovX = double(view.fovX);
            fovY = double(view.fovY);
            displayedSize = size(get(imageHandles(1), 'CData'));
            pixelX = fovX ./ max(displayedSize(2), 1);
            pixelY = fovY ./ max(displayedSize(1), 1);
            set(imageHandles, ...
                'XData', [-fovX/2 + pixelX/2, fovX/2 - pixelX/2], ...
                'YData', [-fovY/2 + pixelY/2, fovY/2 - pixelY/2]);
            app.axis_image.XLim = [-fovX/2, fovX/2];
            app.axis_image.YLim = [-fovY/2, fovY/2];
            axis(app.axis_image, 'image');
        end
        
        function autoWindowFromImage(app, img)
    x = double(img(:));
    x = x(isfinite(x));

    if isempty(x)
        return;
    end

    pLo = prctile(x, 1);
    pHi = prctile(x, 99);

    if pHi <= pLo
        pLo = min(x);
        pHi = max(x);
        if pHi <= pLo
            pHi = pLo + 1;
        end
    end

    WL = (pLo + pHi)/2;
    WW = max(pHi - pLo, eps);

    % Set slider limits wide enough that you can “push” past auto-window
    app.WindowLevelSlider.Limits = [pLo - 0.5*WW, pHi + 0.5*WW];  % WL range
    app.WindowWidthSlider.Limits   = [WW*0.01, WW*5];               % WW range

    app.WindowLevelSlider.Value = WL;
    app.WindowWidthSlider.Value   = WW;
end 
        
function update3DFFTDropDown(app)

    obj = app.scan_data{app.currentscan};

    nSlicesAcquired = obj.slices;

    % Useful through-plane zero-fill choices: acquired size, next power of
    % two, and 2x acquired size. The previous +1/+2 choices added almost no
    % interpolation while still changing odd/even centring.
    fft3dVals = unique([nSlicesAcquired, ...
        2^nextpow2(nSlicesAcquired), 2*nSlicesAcquired], 'stable');

    app.FFTDropDown_3D.Items = cellstr(string(fft3dVals));
    app.FFTDropDown_3D.ItemsData = fft3dVals;

    % Initialise obj.fft_size_3d if missing/empty
    needsDefault = false;

    try
        if ~isprop(obj, 'fft_size_3d') || isempty(obj.fft_size_3d)
            needsDefault = true;
        end
    catch
        try
            if ~isprop(obj, 'fft_size_3d') || isempty(obj.fft_size_3d)
                needsDefault = true;
            end
        catch
            needsDefault = true;
        end
    end

    if needsDefault
        obj.fft_size_3d = nSlicesAcquired;
        app.scan_data{app.currentscan} = obj;
    end

    % If the stored value is not one of the default options, include it
    currentVal = obj.fft_size_3d;

    if ~ismember(currentVal, fft3dVals)
        fft3dVals = unique([fft3dVals, currentVal], 'stable');
        app.FFTDropDown_3D.Items = cellstr(string(fft3dVals));
        app.FFTDropDown_3D.ItemsData = fft3dVals;
    end

    app.FFTDropDown_3D.Value = currentVal;

    % Only enable for 3D scans
    if obj.TwoDimensional == 1 || nSlicesAcquired <= 1
        app.FFTDropDown_3D.Enable = 'off';
        app.Sizeof3DFFTLabel.Enable = 'off';
    else
        app.FFTDropDown_3D.Enable = 'on';
        app.Sizeof3DFFTLabel.Enable = 'on';
    end

end

function updateSliceDropDownFromImages(app)

    [displayScan, displayIndex, isLocalizer] = getDisplayScan(app);
    planeName = currentViewOrientation(app);
    geomT = currentDisplayTransform(app);
    view = ffc_anatomical_plane_geometry(displayScan, planeName, 1, geomT);
    nSlicesDisplay = view.nSlices;

    oldSlice = 1;
    if ~isempty(app.SliceDropDown.Value)
        oldSlice = app.SliceDropDown.Value;
    end

    app.SliceDropDown.Items = cellstr(string(1:nSlicesDisplay));
    app.SliceDropDown.ItemsData = 1:nSlicesDisplay;
    selectedSlice = min(oldSlice, nSlicesDisplay);
    svsIndex = findSVSOverlayScan(app);
    if svsIndex > 0 && scanIsImaging(app, displayScan)
        svsGeometry = ffc_scan_spatial_geometry( ...
            app.scan_data{svsIndex});
        selectedSlice = ffc_nearest_anatomical_slice(displayScan, ...
            planeName, svsGeometry.centerPatientMm, geomT);
    end
    app.SliceDropDown.Value = selectedSlice;

    if nSlicesDisplay <= 1
        app.SliceDropDown.Enable = 'off';
    else
        app.SliceDropDown.Enable = 'on';
    end

end

        function drawimages(app)
            if isempty(app.scan_data) || isempty(app.currentscan)
                return;
            end

            currentScan = app.scan_data{app.currentscan};
            [displayScan, displayIndex, isLocalizer] = getDisplayScan(app);
            if isempty(displayScan.magimage), return; end
            currentIsSVS = scanIsSVS(app, currentScan);
            noLocalizer = currentIsSVS && displayIndex == app.currentscan;

            if noLocalizer
                cla(app.axis_image);
                cla(app.axis_kspace);
                text(app.axis_image, 0.5, 0.5, ...
                    {'SVS relaxometry data loaded', ...
                     'Load an imaging acquisition to display the voxel'}, ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center');
                view = ffc_anatomical_plane_geometry(currentScan, ...
                    acquiredPlaneName(app, currentScan), 1, ...
                    currentDisplayTransform(app));
                drawOrientationCompass(app, view);
                drawnow;
                return;
            end

            dispField = app.EvolutionFieldDropDown.Value;
            dispTime = app.EvolutionTimeDropDown.Value;
            if isLocalizer
                displayField = 1;
                displayTime = 1;
            else
                displayField = min(dispField, size(displayScan.magimage, 5));
                displayTime = min(dispTime, size(displayScan.magimage, 4));
            end

            planeName = currentViewOrientation(app);
            sliceIndex = app.SliceDropDown.Value;
            geomT = currentDisplayTransform(app);
            imageType = string(app.ImageTypeDropDown.Value);
            if currentIsSVS
                imageType = "Magnitude";
            end

            if imageType == "Relaxometry Map" && ...
                    (isempty(currentScan.T1Maps) || isempty(currentScan.R1Maps))
                imageType = "Magnitude";
                app.ImageTypeDropDown.Value = 'Magnitude';
            end
            switch imageType
                case "Magnitude"
                    volume = displayScan.magimage(:, :, :, ...
                        displayTime, displayField);
                case "Phase"
                    volume = displayScan.phaseimage(:, :, :, ...
                        displayTime, displayField);
                otherwise
                    currentScan.R1T1 = app.DisplayingasT1Button.Value;
                    app.scan_data{app.currentscan} = currentScan;
                    if currentScan.R1T1 == 0
                        volume = currentScan.T1Maps(:, :, :, dispField);
                    else
                        volume = currentScan.R1Maps(:, :, :, dispField);
                    end
                    if app.LogarithmicMapsButton.Value
                        volume = log(volume);
                    end
            end

            [displayedImage, view] = ffc_extract_anatomical_plane( ...
                volume, displayScan, planeName, sliceIndex, geomT);
            if imageType == "Magnitude"
                colormap(app.axis_image,gray);
                imageClim = windowLimits(app);
            elseif imageType == "Phase"
                colormap(app.axis_image,gray);
                imageClim = [-pi pi];
            else
                colormap(app.axis_image,jet);
                imageClim = [app.ColourAxisLowerLimitEditField.Value, ...
                    app.ColourAxisUpperLimitEditField.Value];
                if imageClim(2) <= imageClim(1), imageClim(2) = imageClim(1)+eps(max(1,abs(imageClim(1)))); end
            end

            % K-space is still generated from the native acquisition product,
            % but it follows the same post-reconstruction anatomical view when
            % a 3D plane is selected.  No display operation touches k-space
            % used by the reconstruction itself.
            kspaceKey = {displayIndex,displayScan.reconstruction_generation, ...
                displayTime,displayField,planeName,sliceIndex,geomT};
            if ~isequaln(app.KspaceKey,kspaceKey) || (isempty(app.KspaceHandle) || ~isgraphics(app.KspaceHandle))
                kspaceVolume = displayScan.magkspace;
                kspaceScan = displayScan; kspaceScan.magimage = kspaceVolume;
                [displayedKspace,~] = ffc_extract_anatomical_plane( ...
                    kspaceVolume(:,:,:,displayTime,displayField),kspaceScan, ...
                    planeName,sliceIndex,geomT);
                displayedKspace = log1p(double(displayedKspace));
                if (isempty(app.KspaceHandle) || ~isgraphics(app.KspaceHandle))
                    app.KspaceHandle = imagesc(app.axis_kspace,displayedKspace);
                    axis(app.axis_kspace,'image'); colormap(app.axis_kspace,gray);
                else
                    set(app.KspaceHandle,'CData',displayedKspace, ...
                        'XData',[1 size(displayedKspace,2)],'YData',[1 size(displayedKspace,1)]);
                    app.axis_kspace.XLim = [0.5 size(displayedKspace,2)+0.5];
                    app.axis_kspace.YLim = [0.5 size(displayedKspace,1)+0.5];
                end
                app.axis_kspace.CLimMode = 'auto';
                app.KspaceKey = kspaceKey;
            end
            planeKey = {displayIndex,displayScan.magnitude_generation,imageType, ...
                displayTime,displayField,planeName,sliceIndex,geomT};
            if ~isequaln(app.ImagePlaneKey,planeKey)
                if ~isempty(app.mask) && isvalid(app.mask), delete(app.mask); end
                app.mask = []; app.ImagePlaneKey = planeKey;
            end
            if (isempty(app.ImageHandle) || ~isgraphics(app.ImageHandle))
                cla(app.axis_image);
                app.ImageHandle = imagesc(app.axis_image,displayedImage);
                app.ImageHandle.HitTest = 'off';
                app.ImageHandle.PickableParts = 'none';
            else
                app.ImageHandle.CData = displayedImage;
            end
            app.axis_image.CLim = imageClim;
            if imageType == "Relaxometry Map", colormap(app.axis_image,jet);
            else, colormap(app.axis_image,gray); end
            applyImageFovGeometry(app,view);
            drawSVSVoxelOverlay(app,view);
            drawCVIContourOverlay(app,view);
            drawOrientationCompass(app,view);
            drawnow limitrate;
        end

        function updateimages(app, force)  %this method reconstructs the images in case of a change of filter or fft size - more time consuming than a simple redraw
            if nargin < 2, force = false; end
            saveViewState(app);
            app.scan_data{app.currentscan}.fft_size = str2double(app.FFTDropDown.Value);
            app.scan_data{app.currentscan}.fft_size_3d = str2double(string(app.FFTDropDown_3D.Value));
            app.scan_data{app.currentscan}.window_size = app.WindowSizeEditField.Value;
            app.scan_data{app.currentscan}.window_function = app.WindowDropDown.Value;
            app.scan_data{app.currentscan}.multichannel_recon = app.MultichannelDropDown.Value;
            app.scan_data{app.currentscan}.combination = app.CombinationDropDown.Value;
            app.scan_data{app.currentscan}.biascorrect = app.BiasFieldDropDown.Value;
            whiteningCheckBox = getNoiseWhiteningCheckBox(app);
            app.scan_data{app.currentscan}.noise_whitening = ...
                app.scan_data{app.currentscan}.n_receivers > 1 && ...
                logical(whiteningCheckBox.Value);
            if force || ~ffc_reconstruction_valid(app.scan_data{app.currentscan})
                startScanTask(app,'reconstruct',struct('force',force));
            else
                drawimages(app);
            end
        end

    function s = mergeStruct(app, s, override)
    if isempty(override) || ~isstruct(override)
        return
    end

    fn = fieldnames(override);
    for k = 1:numel(fn)
        if ~isempty(override.(fn{k}))
            s.(fn{k}) = override.(fn{k});
        end
    end
end

function redrawCurrentMagImage(app, doAutoWindow)

    if nargin < 2
        doAutoWindow = true;
    end

    if doAutoWindow
        [img0, ~, ~] = currentMagnitudePlane(app);
        if ~isempty(img0)
            autoWindowFromImage(app, img0);
        end
    end

    drawimages(app);

end

function updateFilterButtonTooltip(app)
    obj = app.scan_data{app.currentscan};
    if isempty(obj.denoise_filter)
        filt = 'None';
    else
        filt = obj.denoise_filter;
    end
    app.FilterSettingsButton.Tooltip = {sprintf('Current denoise filter: %s', filt)};
end
    
function openFilterSubGui(app)

    closeFilterDialogs(app);
    targetScanIndex = app.currentscan;
    targetDataset = app.DatasetGeneration;
    obj = app.scan_data{targetScanIndex};
    targetGeneration = obj.magnitude_generation;
    baseMag = obj.filter_source_magimage;
    if isempty(baseMag), baseMag = obj.magimage; end
    previewRunner = FFCJobRunner();
    previewTimer = [];
    cachedKey = []; cachedFiltered = [];
    pendingKey = [];
    origHandle = []; filtHandle = [];

    % Start with the same windowing as the main GUI
    previewWL = app.WindowLevelSlider.Value;
    previewWW = app.WindowWidthSlider.Value;

    % Store last calculated preview so WL/WW changes do not rerun the filter
    lastPreviewImg = [];
    lastFilteredPreview = [];

    % Use stored params if they exist, otherwise sensible defaults.
    % The filter itself is still chosen inside this sub-GUI.
    if isstruct(obj.denoise_params)
        currentParams = mergeStruct(app, ffc_default_denoise_params('LLR+Wavelet'), obj.denoise_params);
    else
        currentParams = ffc_default_denoise_params('LLR+Wavelet');
    end

    filterItems = {'None', ...
                   'Wavelet', ...
                   'LLR', ...
                   'LLR+Wavelet', ...
                   'Total Variation', ...
                   'Generalised Total Variation', ...
                   'LLR+Wavelet+TGV'};

    defaultFilterInSubGui = obj.denoise_filter;
    if isempty(defaultFilterInSubGui) || ~ismember(defaultFilterInSubGui,filterItems)
        defaultFilterInSubGui = 'None';
    end

    sz = size(baseMag);
    if numel(sz) < 5
        sz = [sz ones(1, 5-numel(sz))];
    end

    nSlices = sz(3);
    nTimes  = sz(4);
    nFields = sz(5);

    dlg = uifigure('Name', 'Filter Settings', ...
        'Position', [100 100 1080 640]);
    app.FilterDialogs = dlg;
    dlg.CloseRequestFcn = @(~,~) closeDialog();
    previewTimer = timer('ExecutionMode','singleShot','StartDelay',0.15, ...
        'TimerFcn',@(~,~) launchPreview());

    axOrig = uiaxes(dlg, 'Position', [30 230 450 370]);
    title(axOrig, 'Original');

    axFilt = uiaxes(dlg, 'Position', [590 230 450 370]);
    title(axFilt, 'Filtered preview');

    % ---------------------------------------------------------------------
    % Filter and stack selection controls
    % ---------------------------------------------------------------------

    uilabel(dlg, 'Position', [30 190 100 22], 'Text', 'Filter type');

    ddType = uidropdown(dlg, ...
        'Position', [130 190 190 22], ...
        'Items', filterItems, ...
        'Value', defaultFilterInSubGui);

    uilabel(dlg, 'Position', [340 190 40 22], 'Text', 'Slice');
    ddSlice = uidropdown(dlg, ...
        'Position', [380 190 60 22], ...
        'Items', cellstr(string(1:nSlices)), ...
        'ItemsData', 1:nSlices, ...
        'Value', min(app.SliceDropDown.Value, nSlices));

    uilabel(dlg, 'Position', [460 190 40 22], 'Text', 'Time');
    ddTime = uidropdown(dlg, ...
        'Position', [500 190 60 22], ...
        'Items', cellstr(string(1:nTimes)), ...
        'ItemsData', 1:nTimes, ...
        'Value', min(app.EvolutionTimeDropDown.Value, nTimes));

    uilabel(dlg, 'Position', [580 190 40 22], 'Text', 'Field');
    ddField = uidropdown(dlg, ...
        'Position', [620 190 60 22], ...
        'Items', cellstr(string(1:nFields)), ...
        'ItemsData', 1:nFields, ...
        'Value', min(app.EvolutionFieldDropDown.Value, nFields));

    % ---------------------------------------------------------------------
    % Windowing controls
    % ---------------------------------------------------------------------

    x = double(baseMag(:,:,ddSlice.Value,ddTime.Value,ddField.Value));
    x = x(:);
    x = x(isfinite(x));

    if isempty(x)
        x = [0; 1];
    end

    pLo = prctile(x, 1);
    pHi = prctile(x, 99);

    if pHi <= pLo
        pLo = min(x);
        pHi = max(x);
        if pHi <= pLo
            pHi = pLo + 1;
        end
    end

    autoWL = (pLo + pHi) / 2;
    autoWW = max(pHi - pLo, eps);

    wlLimits = [pLo - 0.5*autoWW, pHi + 0.5*autoWW];
    wwLimits = [autoWW*0.01, autoWW*5];

    previewWL = min(max(previewWL, wlLimits(1)), wlLimits(2));
    previewWW = min(max(previewWW, wwLimits(1)), wwLimits(2));

    uilabel(dlg, 'Position', [710 195 90 22], 'Text', 'Window level');
    slWL = uislider(dlg, ...
        'Position', [800 205 230 3], ...
        'Limits', wlLimits, ...
        'Value', previewWL);

    uilabel(dlg, 'Position', [710 160 90 22], 'Text', 'Window width');
    slWW = uislider(dlg, ...
        'Position', [800 170 230 3], ...
        'Limits', wwLimits, ...
        'Value', previewWW);

    uilabel(dlg,'Position',[710 122 90 22],'Text','LLR grouping');
    grouping = 'all';
    if isfield(currentParams,'stack_grouping'), grouping = currentParams.stack_grouping; end
    if ~ismember(grouping,{'all','per-slice'}), grouping = 'all'; end
    ddGrouping = uidropdown(dlg,'Position',[800 122 230 22], ...
        'Items',{'All slices (legacy)','Per slice'},'ItemsData',{'all','per-slice'}, ...
        'Value',grouping,'Tooltip','Preview and Apply use the same grouping. Per slice changes the LLR model.');
    ddGrouping.ValueChangedFcn = @(~,~) updatePreview();

    % ---------------------------------------------------------------------
    % Denoising parameter controls
    % ---------------------------------------------------------------------

    uilabel(dlg, 'Position', [30 145 80 22], 'Text', 'Patch');
    spPatch = uispinner(dlg, ...
        'Position', [90 145 60 22], ...
        'Limits', [3 31], ...
        'Step', 1, ...
        'Value', currentParams.patch);

    uilabel(dlg, 'Position', [170 145 80 22], 'Text', 'Step');
    spStep = uispinner(dlg, ...
        'Position', [220 145 60 22], ...
        'Limits', [1 16], ...
        'Step', 1, ...
        'Value', currentParams.step);

    uilabel(dlg, 'Position', [300 145 80 22], 'Text', 'LLR iter');
    spNIter = uispinner(dlg, ...
        'Position', [360 145 60 22], ...
        'Limits', [1 10], ...
        'Step', 1, ...
        'Value', currentParams.niter);

    uilabel(dlg, 'Position', [440 145 80 22], 'Text', 'Tau');
    edTau = uieditfield(dlg, 'numeric', ...
        'Position', [480 145 70 22], ...
        'Value', currentParams.tau);

    uilabel(dlg, 'Position', [570 145 90 22], 'Text', 'Wavelet lvl');
    spWav = uispinner(dlg, ...
        'Position', [650 145 60 22], ...
        'Limits', [1 6], ...
        'Step', 1, ...
        'Value', currentParams.wavelet_level);

    uilabel(dlg, 'Position', [30 110 80 22], 'Text', 'TGV a0');
    edA0 = uieditfield(dlg, 'numeric', ...
        'Position', [90 110 70 22], ...
        'Value', currentParams.tgv_alpha0);

    uilabel(dlg, 'Position', [180 110 80 22], 'Text', 'TGV a1');
    edA1 = uieditfield(dlg, 'numeric', ...
        'Position', [240 110 70 22], ...
        'Value', currentParams.tgv_alpha1);

    uilabel(dlg, 'Position', [330 110 90 22], 'Text', 'TGV lambda');
    edLam = uieditfield(dlg, 'numeric', ...
        'Position', [410 110 70 22], ...
        'Value', currentParams.tgv_lambda);

    uilabel(dlg, 'Position', [500 110 90 22], 'Text', 'TGV iter');
    spTGVIter = uispinner(dlg, ...
        'Position', [570 110 70 22], ...
        'Limits', [10 1000], ...
        'Step', 10, ...
        'Value', currentParams.tgv_niter);

    % ---------------------------------------------------------------------
    % Buttons
    % ---------------------------------------------------------------------
    btnAutoTGV = uibutton(dlg, 'push', ...
    'Position', [800 90 150 32], ...
    'Text', 'Auto TGV');


    btnPreview = uibutton(dlg, 'push', ...
        'Position', [30 55 130 32], ...
        'Text', 'Preview');

    btnApply = uibutton(dlg, 'push', ...
        'Position', [180 55 180 32], ...
        'Text', 'Apply to current set');

    btnPermanent = uibutton(dlg, 'push', ...
        'Position', [380 55 230 32], ...
        'Text', 'Reset filter');

    btnResetParams = uibutton(dlg, 'push', ...
        'Position', [630 55 150 32], ...
        'Text', 'Reset parameters');

    btnMatchMainWindow = uibutton(dlg, 'push', ...
        'Position', [800 55 150 32], ...
        'Text', 'Match main window');

    btnClose = uibutton(dlg, 'push', ...
        'Position', [960 55 80 32], ...
        'Text', 'Close');

    statusLabel = uilabel(dlg, ...
        'Position', [30 20 1000 22], ...
        'Text', 'Choose filter and press Preview');

    % ---------------------------------------------------------------------
    % Callbacks
    % ---------------------------------------------------------------------

    ddType.ValueChangedFcn   = @(~,~) filterTypeChanged();
    ddSlice.ValueChangedFcn  = @(~,~) updatePreview();
    ddTime.ValueChangedFcn   = @(~,~) updatePreview();
    ddField.ValueChangedFcn  = @(~,~) updatePreview();

    btnPreview.ButtonPushedFcn         = @(~,~) updatePreview();
    btnApply.ButtonPushedFcn           = @(~,~) applyToCurrentSet();
    btnPermanent.ButtonPushedFcn       = @(~,~) permanentlyApply();
    btnResetParams.ButtonPushedFcn     = @(~,~) resetFilterParameters();
    btnMatchMainWindow.ButtonPushedFcn = @(~,~) matchMainWindow();
    btnClose.ButtonPushedFcn           = @(~,~) closeDialog();
    btnAutoTGV.ButtonPushedFcn = @(~,~) autoTGVParams();

    % Windowing should not rerun the filter
    slWL.ValueChangedFcn = @(~,~) committedWindow();
    slWW.ValueChangedFcn = @(~,~) committedWindow();

    % More responsive while dragging, if your MATLAB version supports it
    slWL.ValueChangingFcn = @(src,event) windowChanging(src, event);
    slWW.ValueChangingFcn = @(src,event) windowChanging(src, event);

    updateControlEnableState();
    updatePreview();

    % =====================================================================
    % Nested functions
    % =====================================================================

    function filterTypeChanged()
        resetFilterParameters();
    end
function autoTGVParams()

    strength = 1.0;  % moderate default

    % Estimate from the same selected-slice stack used for preview
    selectedSlice = ddSlice.Value;
    estStack = baseMag(:,:,selectedSlice,:,:);

    p = ffc_estimate_tgv_params(estStack, strength);

    edA0.Value = p.tgv_alpha0;
    edA1.Value = p.tgv_alpha1;
    edLam.Value = p.tgv_lambda;
    spTGVIter.Value = p.tgv_niter;

    statusLabel.Text = sprintf( ...
        'Auto TGV: sigmaRel %.4g, a0 %.4g, a1 %.4g, lambda %.4g', ...
        p.sigmaRel, p.tgv_alpha0, p.tgv_alpha1, p.tgv_lambda);

    updatePreview();

end
    function [filterType, params, previewStack, previewIndex, key] = readControls()

        filterType = ddType.Value;

        params = struct();

        params.patch = spPatch.Value;
        if mod(params.patch, 2) == 0
            params.patch = params.patch + 1;
            spPatch.Value = params.patch;
        end

        params.step = spStep.Value;
        params.niter = spNIter.Value;
        params.tau = edTau.Value;
        params.sigma = [];
        params.stack_grouping = ddGrouping.Value;
        params.use_parallel = false;
        params.wavelet_level = spWav.Value;

        params.tgv_alpha0 = edA0.Value;
        params.tgv_alpha1 = edA1.Value;
        params.tgv_lambda = edLam.Value;
        params.tgv_niter = spTGVIter.Value;

params.tv_lambda = edLam.Value;
params.tv_niter = spTGVIter.Value;

        selectedSlice = ddSlice.Value;
        selectedTime  = ddTime.Value;
        selectedField = ddField.Value;

        if contains(lower(filterType),'llr')
            if strcmp(params.stack_grouping,'all')
                previewStack = baseMag;
                previewIndex = sub2ind([nSlices nTimes nFields],selectedSlice,selectedTime,selectedField);
                scope = [0 0 0];
            else
                previewStack = baseMag(:,:,selectedSlice,:,:);
                previewIndex = sub2ind([nTimes nFields],selectedTime,selectedField);
                scope = [selectedSlice 0 0];
            end
        else
            previewStack = baseMag(:,:,selectedSlice,selectedTime,selectedField);
            previewIndex = 1;
            scope = [selectedSlice selectedTime selectedField];
        end
        key = {filterType,params,scope};
    end

    function updatePreview()
        if ~isgraphics(dlg), return; end
        updateControlEnableState();
        [~,~,previewStack,previewIndex,key] = readControls();
        originalStack = reshape(previewStack,sz(1),sz(2),[]);
        lastPreviewImg = originalStack(:,:,previewIndex);
        if isequaln(key,cachedKey)
            stop(previewTimer);
            previewRunner.cancel(); pendingKey = [];
            lastFilteredPreview = cachedFiltered(:,:,previewIndex);
            updatePreviewWindowOnly();
            statusLabel.Text = 'Preview ready (cached)';
        elseif ~isequaln(key,pendingKey)
            previewRunner.cancel(); pendingKey = key;
            lastFilteredPreview = []; % never label an older parameter result as current
            if ~isempty(filtHandle) && isgraphics(filtHandle), filtHandle.Visible = 'off'; end
            stop(previewTimer); start(previewTimer);
            statusLabel.Text = 'Preview queued...';
        end
    end

    function launchPreview()
        if ~isgraphics(dlg), return; end
        [filterType,params,previewStack,~,key] = readControls();
        pendingKey = key;
        statusLabel.Text = 'Calculating preview...';
        previewRunner.start('preview',previewStack,struct('filter',filterType,'params',params), ...
            @(result) finishPreview(result,key),@(err) previewFailed(err));
    end

    function finishPreview(result,key)
        if ~isgraphics(dlg), return; end
        [~,~,~,~,currentKey] = readControls();
        pendingKey = [];
        if ~isequaln(key,currentKey), updatePreview(); return; end
        cachedFiltered = reshape(result,sz(1),sz(2),[]);
        cachedKey = key;
        updatePreview();
    end

    function previewFailed(err)
        if ~isgraphics(dlg), return; end
        pendingKey = [];
        statusLabel.Text = ['Preview error: ' err.message];
    end

    function updatePreviewWindowOnly()
        if isempty(lastPreviewImg), return; end
        if isempty(origHandle) || ~isgraphics(origHandle)
            origHandle = imagesc(axOrig,lastPreviewImg);
            axis(axOrig,'image'); colormap(axOrig,gray);
            axOrig.XTick = []; axOrig.YTick = [];
        else
            origHandle.CData = lastPreviewImg;
        end
        if ~isempty(lastFilteredPreview)
            if isempty(filtHandle) || ~isgraphics(filtHandle)
                filtHandle = imagesc(axFilt,lastFilteredPreview);
                axis(axFilt,'image'); colormap(axFilt,gray);
                axFilt.XTick = []; axFilt.YTick = [];
            else
                filtHandle.CData = lastFilteredPreview; filtHandle.Visible = 'on';
            end
        end
        width = max(previewWW,eps(max(1,abs(previewWL))));
        limits = previewWL + [-0.5 0.5]*width;
        axOrig.CLim = limits; axFilt.CLim = limits;
        title(axOrig,sprintf('Source slice %d',ddSlice.Value));
        title(axFilt,sprintf('Preview: %s',ddType.Value));
        drawnow limitrate;
    end

    function committedWindow()
        previewWL = slWL.Value; previewWW = slWW.Value;
        updatePreviewWindowOnly();
    end

    function windowChanging(src,event)
        src.Value = event.Value;
        if src == slWL, previewWL = event.Value;
        else, previewWW = event.Value; end
        updatePreviewWindowOnly();
    end

    function resetFilterParameters()

        filterType = ddType.Value;
        p = ffc_default_denoise_params(filterType);

        spPatch.Value = p.patch;
        spStep.Value = p.step;
        spNIter.Value = p.niter;
        edTau.Value = p.tau;
        spWav.Value = p.wavelet_level;

        if isfield(p, 'tgv_alpha0')
            edA0.Value = p.tgv_alpha0;
        end

        if isfield(p, 'tgv_alpha1')
            edA1.Value = p.tgv_alpha1;
        end

        if isfield(p, 'tgv_lambda')
            edLam.Value = p.tgv_lambda;
        end

        if isfield(p, 'tgv_niter')
            spTGVIter.Value = p.tgv_niter;
        end

        updateControlEnableState();
        updatePreview();

    end

    function matchMainWindow()

        previewWL = app.WindowLevelSlider.Value;
        previewWW = app.WindowWidthSlider.Value;

        previewWL = min(max(previewWL, slWL.Limits(1)), slWL.Limits(2));
        previewWW = min(max(previewWW, slWW.Limits(1)), slWW.Limits(2));

        slWL.Value = previewWL;
        slWW.Value = previewWW;

        updatePreviewWindowOnly();

    end

    function applyToCurrentSet()
        if ~targetIsCurrent(), return; end
        [filterType,params] = readControls();
        closeDialog();
        startScanTask(app,'filter',struct('filter',filterType,'params',params));
    end

    function permanentlyApply()
        % Reset to the same unfiltered/registered source; no cumulative pass.
        if ~targetIsCurrent(), return; end
        closeDialog();
        startScanTask(app,'filter',struct('filter','None', ...
            'params',ffc_default_denoise_params('None')));
    end

    function valid = targetIsCurrent()
        valid = targetDataset == app.DatasetGeneration && ...
            targetScanIndex == app.currentscan && ...
            targetGeneration == app.scan_data{targetScanIndex}.magnitude_generation;
        if ~valid
            statusLabel.Text = 'The source changed. Close and reopen Filter Settings.';
        end
    end

    function closeDialog()
        if ~isempty(previewTimer) && isvalid(previewTimer)
            stop(previewTimer); delete(previewTimer); previewTimer = [];
        end
        if isvalid(previewRunner), delete(previewRunner); end
        if isgraphics(dlg), delete(dlg); end
    end

    function updateControlEnableState()

        filterType = ddType.Value;

        usesLLR = ismember(filterType, ...
            {'LLR','LLR+Wavelet','LLR+Wavelet+TGV'});

        usesWavelet = ismember(filterType, ...
            {'Wavelet','LLR+Wavelet','LLR+Wavelet+TGV'});

        usesTGV = ismember(filterType, ...
            {'Generalised Total Variation','LLR+Wavelet+TGV'});

        usesTV = ismember(filterType, ...
            {'Total Variation'});

        llrState = matlab.lang.OnOffSwitchState(usesLLR);
        wavState = matlab.lang.OnOffSwitchState(usesWavelet);
        tgvState = matlab.lang.OnOffSwitchState(usesTGV || usesTV);

        ddGrouping.Enable = llrState;
        spPatch.Enable = llrState;
        spStep.Enable = llrState;
        spNIter.Enable = llrState;
        edTau.Enable = llrState;

        spWav.Enable = wavState;

        edA0.Enable = tgvState;
        edA1.Enable = tgvState;
        edLam.Enable = tgvState;
        spTGVIter.Enable = tgvState;

    end

end   


    end
    
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.JobRunner = FFCJobRunner();
            app.UIFigure.Name = 'FFC Processing - Performance Edition + CVI contours TEST 2';

            cviMenu = uimenu(app.UIFigure,'Text','CVI Contours');
            app.cviMenuHandles.import = uimenu(cviMenu,'Text','Import workspace...', ...
                'MenuSelectedFcn',@(src,event) importCVIContours(app));
            app.cviMenuHandles.show = uimenu(cviMenu,'Text','Show contours', ...
                'Checked','off','MenuSelectedFcn',@(src,event) toggleCVIContours(app));
            app.cviMenuHandles.registration = uimenu(cviMenu, ...
                'Text','Registration / display...', ...
                'MenuSelectedFcn',@(src,event) openCVIRegistrationDialog(app));
            app.cviMenuHandles.clear = uimenu(cviMenu,'Text','Clear contours', ...
                'Separator','on','MenuSelectedFcn',@(src,event) clearCVIContours(app));

            % Keep the contour workflow visible in the main window as well as
            % in the menu.  This avoids relying on platform-specific menu-bar
            % presentation and makes it obvious which build is running.
            app.ScanList.Position = [19 43 260 228];
            app.cviMenuHandles.importButton = uibutton(app.ScanListPanel,'push', ...
                'Text','Import CVI contours', ...
                'Position',[19 11 124 24], ...
                'Tooltip',{'Import ventricular contours from a cvi42 workspace'}, ...
                'ButtonPushedFcn',@(src,event) importCVIContours(app));
            app.cviMenuHandles.registrationButton = uibutton(app.ScanListPanel,'push', ...
                'Text','CVI registration', ...
                'Position',[153 11 126 24], ...
                'Tooltip',{'Open the CVI-to-FCI rigid registration controls'}, ...
                'ButtonPushedFcn',@(src,event) openCVIRegistrationDialog(app));

            app.WindowLevelSlider.ValueChangingFcn = @(src,event) mainWindowChanging(app,src,event);
            app.WindowWidthSlider.ValueChangingFcn = @(src,event) mainWindowChanging(app,src,event);
            % This control is created at run time so the existing App Designer
            % component model remains backwards compatible.
            whiteningCheckBox = uicheckbox(app.panel_recon);
            whiteningCheckBox.Text = 'Noise whitening';
            whiteningCheckBox.Value = false;
            whiteningCheckBox.Enable = 'off';
            whiteningCheckBox.Position = [11 68 145 22];
            whiteningCheckBox.Tooltip = {'Noise whitening requires multiple receive channels.'};
            whiteningCheckBox.ValueChangedFcn = ...
                @(source, event) noiseWhiteningValueChanged(app, source);
            setappdata(app.panel_recon, 'NoiseWhiteningCheckBox', whiteningCheckBox);

            orientationLabel = uilabel(app.CurrentImagePanel);
            orientationLabel.Text = 'View';
            orientationLabel.HorizontalAlignment = 'right';
            orientationLabel.Position = [230 15 48 22];
            orientationLabel.Enable = 'off';
            orientationDropDown = uidropdown(app.CurrentImagePanel);
            orientationDropDown.Items = {'Axial', 'Coronal', 'Sagittal'};
            orientationDropDown.Value = 'Axial';
            orientationDropDown.Position = [293 15 100 22];
            orientationDropDown.Enable = 'off';
            orientationDropDown.ValueChangedFcn = ...
                @(source, event) viewOrientationValueChanged(app, source);
            setappdata(app.CurrentImagePanel, ...
                'ViewOrientationDropDown', orientationDropDown);
            setappdata(app.CurrentImagePanel, ...
                'ViewOrientationLabel', orientationLabel);
        end

        % Button pushed function: LoadDataSetButton
        function LoadDataSetButtonPushed(app, event)
            cancelScanTask(app);
            closeFilterDialogs(app);
            file = open_mrifile;
            if isempty(file)
                return; % cancellation must not replace an existing dataset
            end
            options = struct( ...
                'backgroundselect', get(app.BackgroundselectionButton, 'Value'), ...
                'reprocess', get(app.ReprocessingButton, 'Value'), ...
                'checkfit', get(app.CheckfitsButton, 'Value'));
            [file, sequencelist, failures] = ffc_preprocess_scans(file, options);
            for n = 1:numel(failures)
                fprintf(2, 'Error while processing sequence number %d:\n%s\n', ...
                    failures(n).sequenceIndex, ...
                    getReport(failures(n).exception, 'extended', 'hyperlinks', 'off'));
            end
            if isempty(file)
                uialert(app.UIFigure, ...
                    ['No sequences could be preprocessed. The current dataset ', ...
                     'was not changed. Full error details are in the Command Window.'], ...
                    'Dataset not loaded');
                return;
            end

            app.DatasetGeneration = app.DatasetGeneration + 1;
            app.ViewStates = {}; app.KspaceKey = []; app.CompassKey = [];
            app.ScanList.Items = sequencelist;
            app.ScanList.ItemsData =1:numel(app.ScanList.Items); %allows us to get index of selected scan
            app.ScanList.Value = 1;
            app.currentscan = 1; %start on scan 1
            app.scan_data = file; %save the scan to the GUI for easy access
            app.svsOverlayScan = [];
            app.svsLocalizerScan = [];
            
            set(findall(app.UIFigure, '-property', 'Enable'), 'Enable', 'on') %enable all the GUI options now that we have a data set to work on. Errors abound if there is no data set loaded and buttons are pushed!

            activeimage(app) %start drawing images
  

            figure(app.UIFigure);
            if ~isempty(failures)
                uialert(app.UIFigure, ...
                    sprintf(['Loaded %d sequence(s); skipped %d failed sequence(s). ', ...
                        'Full error details are in the Command Window.'], ...
                        numel(file), numel(failures)), ...
                    'Some sequences were not loaded', 'Icon', 'warning');
            end
            
        end

        % Value changed function: WindowLevelSlider
        function WindowLevelSliderValueChanged(app, event)
            updateWindowOnly(app)
        end

        % Value changed function: WindowWidthSlider
        function WindowWidthSliderValueChanged(app, event)
          updateWindowOnly(app)
        end

        % Callback function
        function GammaSliderValueChanged(app, event)
        drawimages(app)
        end

        % Value changed function: ScanList
        function ScanListValueChanged(app, event)
            saveViewState(app);
            cancelScanTask(app);
            closeFilterDialogs(app);
            activeimage(app)
            
        end

        % Value changed function: SliceDropDown
        function SliceDropDownValueChanged(app, event)
            drawimages(app)
        end

        function viewOrientationValueChanged(app, source)
            if isempty(app.scan_data) || isempty(app.currentscan)
                return;
            end
            currentScan = app.scan_data{app.currentscan};
            currentScan.view_orientation = char(string(source.Value));
            app.scan_data{app.currentscan} = currentScan;
            updateSliceDropDownFromImages(app);
            [image, ~, ~] = currentMagnitudePlane(app);
            if ~isempty(image)
                autoWindowFromImage(app, image);
            end
            drawimages(app);
        end

        % Value changed function: FFTDropDown
        function FFTDropDownValueChanged(app, event)
            updateimages(app)
        end

        % Value changed function: EvolutionTimeDropDown
        function EvolutionTimeDropDownValueChanged(app, event)
            drawimages(app)
        end

        % Value changed function: EvolutionFieldDropDown
        function EvolutionFieldDropDownValueChanged(app, event)
            disp_field = app.EvolutionFieldDropDown.Value;
            tevol_list = sprintfc('%g',round(app.scan_data{app.currentscan}.timepoints(disp_field,:),2)); %update the evolution times
            app.EvolutionTimeDropDown.Items = tevol_list;
            app.EvolutionTimeDropDown.ItemsData = 1:numel(app.EvolutionTimeDropDown.Items);
            drawimages(app)
        end

        % Value changed function: ImageTypeDropDown
        function ImageTypeDropDownValueChanged(app, event)
            drawimages(app)
        end

        % Value changed function: WindowDropDown
        function WindowDropDownValueChanged(app, event)
            updateimages(app)
        end

        % Value changed function: WindowSizeEditField
        function WindowSizeEditFieldValueChanged(app, event)
            updateimages(app)
        end

        % Button pushed function: GenerateROIDispersionCurveButton
        function GenerateROIDispersionCurveButtonPushed(app, event)

    targetIndex = app.currentscan; targetDataset = app.DatasetGeneration;
    targetGeneration = app.scan_data{targetIndex}.magnitude_generation;
    currentScan = app.scan_data{app.currentscan};
    currentScan.checkfit = app.CheckfitsButton.Value;
    currentScan.R1T1 = app.DisplayingasT1Button.Value;

    if scanIsSVS(app, currentScan)
        % There is no spatial ROI to draw for SVS. buildSVSSignal has
        % already reduced dimension 1 using mean(abs(rawdata),1).
        currentScan.mask = true(1, 1);
        [currentScan, ~] = T1dispersion(currentScan, 1);
        app.scan_data{app.currentscan} = currentScan;
        return;
    end

    currentView = ffc_anatomical_plane_geometry(currentScan, ...
        currentViewOrientation(app), app.SliceDropDown.Value, ...
        currentDisplayTransform(app));
    if currentView.normalDimension ~= 3
        acquiredPlane = acquiredPlaneName(app, currentScan);
        uialert(app.UIFigure, ...
            sprintf(['ROI dispersion uses the acquired image planes. ', ...
                'Switch View to %s before drawing the ROI.'], acquiredPlane), ...
            'Select acquired plane');
        return;
    end

    % Delete previous ROI safely
    try
        if ~isempty(app.mask) && isvalid(app.mask)
            delete(app.mask);
        end
    catch
    end

    % --- Temporarily disable custom mouse handlers so drawfreehand can work ---
    f = app.UIFigure;

    oldDown  = f.WindowButtonDownFcn;
    oldMove  = f.WindowButtonMotionFcn;
    oldUp    = f.WindowButtonUpFcn;

    % Scroll wheel is optional depending on what you added
    oldWheel = [];
    if isprop(f,'WindowScrollWheelFcn')
        oldWheel = f.WindowScrollWheelFcn;
    end

    f.WindowButtonDownFcn = [];
    f.WindowButtonMotionFcn = [];
    f.WindowButtonUpFcn = [];
    if isprop(f,'WindowScrollWheelFcn')
        f.WindowScrollWheelFcn = [];
    end

    restoreMouse = onCleanup(@() restoreMouseHandlers(app,f,oldDown,oldMove,oldUp,oldWheel));
    % Give focus to the figure (helps interaction)
    figure(f);

    % --- Draw ROI (blocking until completed/double-click) ---
    app.mask = drawfreehand(app.axis_image);

    clear restoreMouse

    % If user cancelled / ROI invalid, bail out cleanly
    if isempty(app.mask) || ~isvalid(app.mask)
        return;
    end

    if targetDataset ~= app.DatasetGeneration || targetIndex ~= app.currentscan || ...
            targetGeneration ~= app.scan_data{targetIndex}.magnitude_generation, return; end
    % Restore the drawn display mask to native reconstruction dimensions.
    currentScan.mask = ffc_restore_native_plane_mask( ...
        createMask(app.mask), currentView);
    [currentScan, ~] = T1dispersion(currentScan, currentView.sliceIndex);
    app.scan_data{app.currentscan} = currentScan;
        end

        % Button pushed function: GenerateRelaxationMapsButton
        function GenerateRelaxationMapsButtonPushed(app, event)
            app.scan_data{app.currentscan}.checkfit = app.CheckfitsButton.Value;
            app.scan_data{app.currentscan}.R1T1 = app.DisplayingasT1Button.Value;
            startScanTask(app,'map',struct());
        end

        % Button pushed function: RotateImageButton
        function RotateImageButtonPushed(app, event)
            Rz = [0 1 0; -1 0 0; 0 0 1];
            [~, displayIndex, ~] = getDisplayScan(app);
            app.scan_data{displayIndex}.geomT = ...
                Rz * app.scan_data{displayIndex}.geomT;
            app.scan_data{displayIndex}.rot_theta = mod( ...
                app.scan_data{displayIndex}.rot_theta + 90, 360);
            drawimages(app)
        end

        % Callback function
        function FilterDropDownValueChanged(app, event)
    
        end

        % Callback function
        function FilterControlEditFieldValueChanged(app, event)
    
            
        end

        % Value changed function: ColourAxisUpperLimitEditField
        function ColourAxisUpperLimitEditFieldValueChanged(app, event)
         %defines the limits of T1 or R1 that are displayed in the
         %relaxation maps
            if app.ColourAxisUpperLimitEditField.Value < app.ColourAxisLowerLimitEditField.Value
                app.ColourAxisUpperLimitEditField.Value = app.ColourAxisLowerLimitEditField.Value + 1; %stops errors when upperbound is less than lowerbound or vice versa
            end
            
            if app.ColourAxisLowerLimitEditField.Value > app.ColourAxisUpperLimitEditField.Value
                app.ColourAxisLowerLimitEditField.Value = app.ColourAxisUpperLimitEditField.Value - 1;
            end
            
            drawimages(app)
        end

        % Value changed function: ColourAxisLowerLimitEditField
        function ColourAxisLowerLimitEditFieldValueChanged(app, event)
            if app.ColourAxisUpperLimitEditField.Value < app.ColourAxisLowerLimitEditField.Value
               app.ColourAxisUpperLimitEditField.Value = app.ColourAxisLowerLimitEditField.Value + 1;
            end
            
            if app.ColourAxisLowerLimitEditField.Value > app.ColourAxisUpperLimitEditField.Value
               app.ColourAxisLowerLimitEditField.Value = app.ColourAxisUpperLimitEditField.Value - 1;
            end
            
            drawimages(app)
        end

        % Value changed function: LogarithmicMapsButton
        function LogarithmicMapsButtonValueChanged(app, event)
            value = app.LogarithmicMapsButton.Value;
            if value
                app.LogarithmicMapsButton.Text = "Log scale Relaxation Maps";
                %try and select sensible limits when log scale selected
                if  app.ColourAxisLowerLimitEditField.Value <1  %log of 0 is inf
                    app.ColourAxisLowerLimitEditField.Value = 1;
                end
                app.ColourAxisLowerLimitEditField.Value =  log(app.ColourAxisLowerLimitEditField.Value);
                app.ColourAxisUpperLimitEditField.Value =  log(app.ColourAxisUpperLimitEditField.Value);
            else
                app.LogarithmicMapsButton.Text = "Linear scale Relaxation Maps";
                app.ColourAxisLowerLimitEditField.Value =  exp(app.ColourAxisLowerLimitEditField.Value);
                app.ColourAxisUpperLimitEditField.Value =  exp(app.ColourAxisUpperLimitEditField.Value);
            end
            drawimages(app)
            
        end

        % Value changed function: DisplayingasT1Button
        function DisplayingasT1ButtonValueChanged(app, event)
            value = app.DisplayingasT1Button.Value;
            if value
                %try and select sensible limits since R1 is in s-1 and T1
                %is in ms
                app.DisplayingasT1Button.Text = "Displaying as R1";
                app.ColourAxisLowerLimitEditField.Value =  app.ColourAxisLowerLimitEditField.Value/100;
                app.ColourAxisUpperLimitEditField.Value =  app.ColourAxisUpperLimitEditField.Value/100;
            else
                app.DisplayingasT1Button.Text = "Displaying as T1";
                app.ColourAxisLowerLimitEditField.Value =  app.ColourAxisLowerLimitEditField.Value*100;
                app.ColourAxisUpperLimitEditField.Value =  app.ColourAxisUpperLimitEditField.Value*100;
            end
            
            drawimages(app)
        end

        % Button pushed function: ExportscantoworkspaceButton
        function ExportscantoworkspaceButtonPushed(app, event)
         file = app.scan_data{app.currentscan};

    try
        rawName = file.displayname;
    catch exception
        rawName = "scan";
    end
        
 

    % Convert to char/string safely
    rawName = string(rawName);

    % Make valid MATLAB variable name
    name = matlab.lang.makeValidName(rawName);

 
    existing = evalin('base','who');
    name = matlab.lang.makeUniqueStrings(name, existing);

    assignin('base', name, file);

    % Optional: tell the user what name it used
    uialert(app.UIFigure, "Exported as: " + name, "Export complete");
    
        end

        % Button pushed function: ExportscanparameterstoworkspaceButton
        function ExportscanparameterstoworkspaceButtonPushed(app, event)
            %extract all the relevant scan parameters from the currently
            %selected file and display them in the command window
            file = app.scan_data{app.currentscan};
            disp('--------------------')
            disp(['Experiment name = ',file.name])
            disp(['Sequence name = ',file.sequence])
            disp(['Matrix size = ',num2str(file.samples),' x ',num2str(file.samples)])
            if file.partialkspace == 1
                disp('Partial Fourier = True')
            else
                disp('Partial Fourier = False')
            end
            disp(['Acquisition Bandwidth = ',num2str(file.bandwidth),' KHz'])
            disp(['FOV = ',num2str(file.fov),' mm'])
            disp(['TE = ',num2str(file.TE),' ms'])
            disp(['Flip Angle = ',num2str(file.FlipAngle),' degrees'])
            disp_slice = app.SliceDropDown.Value;
            orientation = file.orientation(:,disp_slice)';
            orientation = num2str(orientation);
            switch orientation
                case '-90   0   0'
                    disp('Orientation = Coronal');
                case '0  0  0'
                    disp('Orientation = Axial');
                otherwise
                    disp(['Orientation = ',orientation]);
            end
            disp(['Slice thickness = ',num2str(file.thk),' mm'])
            disp(['Inplane resolution = = ',num2str(file.resolution_inplane),' mm'])
            disp('Evolution Fields (mT): ')
            disp(num2str(file.fieldpoints,3))
            disp('Evolution Times (ms): ')
            disp(num2str(file.timepoints,3))
            
            
            disp('--------------------')
        end

        % Button pushed function: kspaceButton
        function kspaceButtonPushed(app, event)
            h = figure;
            %get the image as it appears in the axis
            A = getimage(app.axis_kspace);
            cmap = colormap(app.axis_kspace);
            colormap(h,cmap)
            imagesc(A); axis off image
        end

        % Button pushed function: ImageButton
        function ImageButtonPushed(app, event)
            if isempty(app.ImageHandle) || ~isgraphics(app.ImageHandle), return; end
            h = figure;
            ax = axes(h);
            imagesc(ax,app.ImageHandle.XData,app.ImageHandle.YData,app.ImageHandle.CData);
            ax.CLim = app.axis_image.CLim;
            colormap(ax,colormap(app.axis_image)); axis(ax,'off','image');
        end

        % Value changed function: IndividualNormalisationButton
        function IndividualNormalisationButtonValueChanged(app, event)
            value = app.IndividualNormalisationButton.Value;
            if value
                app.IndividualNormalisationButton.Text = "Relative Normalisation";
            else
                app.IndividualNormalisationButton.Text = "Individual Normalisation";
            end
            drawimages(app)
        end

        % Button pushed function: FlipHorizontallyButton
        function FlipHorizontallyButtonPushed(app, event)
            Sx = [-1 0 0; 0 1 0; 0 0 1];
            [~, displayIndex, ~] = getDisplayScan(app);
            app.scan_data{displayIndex}.geomT = ...
                Sx * app.scan_data{displayIndex}.geomT;
            drawimages(app)
           
        end

        % Button pushed function: FlipVerticallyButton
        function FlipVerticallyButtonPushed(app, event)
            Sy = [1 0 0; 0 -1 0; 0 0 1];
            [~, displayIndex, ~] = getDisplayScan(app);
            app.scan_data{displayIndex}.geomT = ...
                Sy * app.scan_data{displayIndex}.geomT;
            drawimages(app)
        end

        % Button pushed function: ExportallimagesButton
        function ExportallimagesButtonPushed(app, event)
            absolute_scale = app.IndividualNormalisationButton.Value;
            ImageType = app.ImageTypeDropDown.Value;
            fieldlist = app.scan_data{app.currentscan}.fieldpoints;
            timepoints =app.scan_data{app.currentscan}.timepoints;
            n_timepoints = app.scan_data{app.currentscan}.n_timepoints;
            n_fieldpoints = app.scan_data{app.currentscan}.n_fieldpoints;
            
            switch ImageType
                case "Relaxometry Map"
                    if app.DisplayingasT1Button.Value
                        imagesforexport = squeeze(app.scan_data{app.currentscan}.R1Maps);
                        tag = 'R1';
                    else
                        imagesforexport = squeeze(app.scan_data{app.currentscan}.T1Maps);
                        tag = 'T1';
                    end
                    
                    for field = 1:n_fieldpoints
                        imagename = char(([num2str(round(fieldlist(field))) ' mT_' tag 'map.tif']));
                        minD = app.ColourAxisLowerLimitEditField.Value;
                        maxD = app.ColourAxisUpperLimitEditField.Value*1;
                        mapped_image = ind2rgb(double(round((1+imagesforexport(:,:,1,field) - minD) * 255 ./ (maxD - minD))), jet(256)); %this creates a tif RGB image with the same colour scaling as in the GUI
                        imwrite(mapped_image, imagename);
                    end
                    
                case "Phase"
                    imagesforexport = squeeze(app.scan_data{app.currentscan}.phaseimage);
                    tag = 'phasemap';
                    for field = 1:n_fieldpoints
                        for times = 1:n_timepoints
                            imagename = char(([num2str(round(fieldlist(field))) ' mT_' num2str(round(timepoints(field,times))),' ms ' tag '.tif']));
                            imwrite(mat2gray(imagesforexport(:,:,times,field)),imagename);
                        end
                    end
                    
                case "Magnitude"
                    imagesforexport = app.scan_data{app.currentscan}.magimage;
                    exportMax = max(imagesforexport(:)) * 0.2;
                    exportMax = max(exportMax,eps('single'));
                   for slice = 1:app.scan_data{app.currentscan}.slices
                    for field = 1:n_fieldpoints
                        for times = 1:n_timepoints
                            imagename = char(([num2str(round(fieldlist(field))) ' mT_' num2str(round(timepoints(field,times).*100)),' ms_slice_' num2str(slice) '.tif']));
                            if absolute_scale
                                imwrite(mat2gray(imagesforexport(:,:,slice,times,field),[0 exportMax]),imagename);
                            else
                                imwrite(mat2gray(imagesforexport(:,:,slice,times,field)),imagename);
                            end
                        end
                    end
                   end
            end
            
        end

        % Button pushed function: ExportstackedimagesButton
        function ExportstackedimagesButtonPushed(app, event)
            absolute_scale = app.IndividualNormalisationButton.Value;
            ImageType = app.ImageTypeDropDown.Value;
            fieldlist = app.scan_data{app.currentscan}.fieldpoints;
            n_timepoints = app.scan_data{app.currentscan}.n_timepoints;
            n_fieldpoints = app.scan_data{app.currentscan}.n_fieldpoints;
            
            switch ImageType
                case "Relaxometry Map"
                    if app.DisplayingasT1Button.Value
                        imagesforexport = squeeze(app.scan_data{app.currentscan}.R1Maps);
                        tag = 'R1';
                    else
                        imagesforexport = squeeze(app.scan_data{app.currentscan}.T1Maps);
                        tag = 'T1';
                    end
                    
                    for field = 1:n_fieldpoints
                        imagename = char((['Stacked ' tag 'map.tif']));
                        minD = app.ColourAxisLowerLimitEditField.Value;
                        maxD = app.ColourAxisUpperLimitEditField.Value;
                        mapped_image = ind2rgb(double(round((1+imagesforexport(:,:,field) - minD) * 255 ./ (maxD - minD))), jet(256)); %this creates a tif RGB image with the same colour scaling as in the GUI
                        imwrite(mapped_image, imagename,'WriteMode','append');
                    end
                    
                case "Phase"
                    imagesforexport = squeeze(app.scan_data{app.currentscan}.phaseimage);
                    for field = 1:n_fieldpoints
                        imagename = char(([num2str(round(fieldlist(field))) ' mT phase stack.tif']));
                        if exist(imagename, 'file')==2
                            delete(imagename);
                        end
                        for times = 1:n_timepoints
                            
                            imwrite(mat2gray(imagesforexport(:,:,times,field)),imagename,'WriteMode','append'); %writemode append adds new images to an existing one, creating a stacked image
                        end
                    end
                    
                case "Magnitude"
                    imagesforexport = app.scan_data{app.currentscan}.magimage;
                    exportMax = max(imagesforexport(:)) * 0.2;
                    exportMax = max(exportMax,eps('single'));
                 %%   for slice = 1:app.scan_data{app.currentscan}.slices
                    for slice = 1:1
                    for field = 1:n_fieldpoints
                        imagename = char(([num2str(round(fieldlist(field))) ' mT stack.tif']));
                        if exist(imagename, 'file')==2
                            delete(imagename);
                        end
                        for times = 1:n_timepoints
                            if absolute_scale
                                imwrite(mat2gray(imagesforexport(:,:,slice,times,field),[0 exportMax]),imagename,'WriteMode','append');
                            else
                                imwrite(mat2gray(imagesforexport(:,:,slice,times,field)),imagename,'WriteMode','append');
                            end
                        end
                    end
                    end
            end
        end

        % Callback function
        function ClusterAnalysisButtonPushed(app, event)
            try FFCAnalysis_V3(app.scan_data{app.currentscan})
            catch ME %#ok<NASGU>
                errordlg('Error launching FFC Analysis. Check that FFCAnalysis is included in the MATLAB path')
            end
        end

        % Button pushed function: RegisterImagesButton
        function RegisterImagesButtonPushed(app, event)
            index = app.currentscan; dataset = app.DatasetGeneration;
            generation = app.scan_data{index}.magnitude_generation;
            try
                masks = ffc_registration_masks(app.scan_data{index}.magimage);
                if dataset ~= app.DatasetGeneration || index ~= app.currentscan || ...
                        generation ~= app.scan_data{index}.magnitude_generation, return; end
                startScanTask(app,'register',struct('masks',masks));
            catch ME
                if ~strcmp(ME.identifier,'ffc_registration_masks:Cancelled')
                    uialert(app.UIFigure,ME.message,'Registration');
                end
            end
        end

        % Button pushed function: SNRButton
        function SNRButtonPushed(app, event)
            %Can be hardcoded once slice known
            disp_slice = app.SliceDropDown.Value;
            disp_field = app.EvolutionFieldDropDown.Value;
            disp_time = app.EvolutionTimeDropDown.Value;

            targetSlice = abs(app.scan_data{app.currentscan}.magimage(:,:,disp_slice,disp_time,disp_field));

            %Dimensions
            rows = app.scan_data{app.currentscan}.views;
            cols = app.scan_data{app.currentscan}.samples;

            %Normalising
            targetSliceNorm = mat2gray(targetSlice);
            maskInitial = imbinarize(targetSliceNorm); %Add this line of code if the mask isnt correct , "adaptive", Sensitivity=0.7);
            se = strel('disk', 4);
            maskFilled = imclose(maskInitial,se);
            
            %Identify and count pixels in connected groups of pixels in the binary mask
            CC = bwconncomp(maskFilled);
            %Applying the numel function to every cell (group of connected pixels) in
            %the CC array
            numPixels = cellfun(@numel, CC.PixelIdxList);
            %Finding the largest region of pixels; the phantom
            [~, idx] = max(numPixels);

            %Create a binary mask with the largest region of pixels
            maskMROI = zeros(rows, cols, 'logical');
            maskMROI(CC.PixelIdxList{idx}) = true;
    
            %Using imerode to shrink the mask to ensure none of the edge pixels are included in the signal average
            seErode = strel('disk', 6);
            maskMROIRefined = imerode(maskMROI, seErode);
            
            %Mean signal intensity
            S = mean(targetSlice(maskMROIRefined));
    
            %Number and size of noise ROI's
            noiseRegions = 4;
            noiseStd = zeros(noiseRegions, 1);
            noiseBoxSize = round(rows * 0.05);
            buffer = 8;
    
            %Creating a table of starting coordinates to place the noise regions, with an 8 pixel buffer from the edges
            cornerCoords = [
                buffer, buffer;
                cols - noiseBoxSize - buffer, buffer;
                cols - noiseBoxSize - buffer, rows - noiseBoxSize - buffer;
                buffer, rows - noiseBoxSize - buffer
            ];

            figure('Name', 'SNR Analysis', 'NumberTitle', 'off');
            overlayImg = labeloverlay(targetSliceNorm, maskMROIRefined, 'Transparency', 0.7, 'Colormap', 'hot');
            imshow(overlayImg);
            title(['SNR Analysis - Slice ', num2str(disp_slice)], 'FontSize', 14);
            axis square
            hold on;
    
            for i = 1:noiseRegions
                xStart = cornerCoords(i, 1);
                yStart = cornerCoords(i, 2);
        
                maskNoise = false(rows, cols);
                maskNoise(yStart : yStart + noiseBoxSize, xStart : xStart + noiseBoxSize) = true;
                noiseStd(i) = std(targetSlice(maskNoise));
        
                rectangle('Position', [xStart, yStart, noiseBoxSize, noiseBoxSize], 'EdgeColor', 'g', 'LineWidth', 2);
            end
            hold off;

            %Calculation of the mean noise which is corrected by dividing by 0.66 which accounts for the Rayleigh distribution of the noise in the magnitude image.
            noiseAvg = mean(noiseStd);
            noiseCorrected = noiseAvg / 0.66; 
            %SNR calculation
            snr = S / noiseCorrected;

            fprintf('--- Signal to Noise Ratio ---\n');
            fprintf('Slice %d SNR: %.3f | Mean Signal: %.4f\n | Mean Noise: %.4f\n', disp_slice, snr, S, noiseAvg);
            
        end

        % Value changed function: MultichannelDropDown
        function MultichannelDropDownValueChanged2(app, event)
          %updateimages(app)
        end

        function noiseWhiteningValueChanged(app, source)
            if isempty(app.scan_data) || isempty(app.currentscan)
                source.Value = false;
                return;
            end

            if app.scan_data{app.currentscan}.n_receivers <= 1
                source.Value = false;
                source.Enable = 'off';
                app.scan_data{app.currentscan}.noise_whitening = false;
                return;
            end

            app.scan_data{app.currentscan}.noise_whitening = logical(source.Value);
            updateimages(app);
        end

        % Button pushed function: RebuildImagesButton
        function RebuildImagesButtonPushed(app, event)
            if app.JobRunner.isRunning()
                cancelScanTask(app);
                if ffc_reconstruction_valid(app.scan_data{app.currentscan}), activeimage(app); end
            else
                updateimages(app,true);
            end 
        end

        % Key press function: UIFigure
        function UIFigureKeyPress(app, event)

    % Ignore keypresses when user is typing into controls
    obj = app.UIFigure.CurrentObject;
    if ~isempty(obj) && isprop(obj,'Type')
        if any(strcmp(obj.Type, {'uieditfield','uitextarea','uidropdown','uispinner'}))
            return;
        end
    end

    switch event.Key

        % -------- Slice navigation --------
        case 'uparrow'
            v = app.SliceDropDown.Value;
            if v > min(app.SliceDropDown.ItemsData)
                app.SliceDropDown.Value = v - 1;
                drawimages(app);
            end

        case 'downarrow'
            v = app.SliceDropDown.Value;
            if v < max(app.SliceDropDown.ItemsData)
                app.SliceDropDown.Value = v + 1;
                drawimages(app);
            end

        % -------- Evolution FIELD navigation --------
        case 'leftarrow'
            v = app.EvolutionTimeDropDown.Value;
            if v > min(app.EvolutionTimeDropDown.ItemsData)
                app.EvolutionTimeDropDown.Value = v - 1;
               drawimages(app);
            end

        case 'rightarrow'
            v = app.EvolutionTimeDropDown.Value;
            if v < max(app.EvolutionTimeDropDown.ItemsData)
                app.EvolutionTimeDropDown.Value = v + 1;
                drawimages(app);
            end
    end


            
        end

        % Button down function: axis_image
        function axis_imageButtonDown(app, event)
            figure(app.UIFigure);   % give UIFigure focus so arrows work
        end

        % Window key press function: UIFigure
        function UIFigureWindowKeyPress(app, event)
            if isempty(app.scan_data), return; end
            UIFigureKeyPress(app,event);
        end

        % Window button down function: UIFigure
        function UIFigureWindowButtonDown(app, event)
        
    % Only respond if click is on the image axes (or an object inside it)
    obj = app.UIFigure.CurrentObject;

    onImageAxes = isequal(obj, app.axis_image);
    if ~onImageAxes
        try
            onImageAxes = isprop(obj,'Parent') && isequal(obj.Parent, app.axis_image);
        catch
        end
    end
    if ~onImageAxes
        return;
    end

    % ---- DOUBLE CLICK = AUTO WINDOW ----
    if strcmp(app.UIFigure.SelectionType, 'open')

        [img, ~, ~] = currentMagnitudePlane(app);
        if ~isempty(img)
            autoWindowFromImage(app, img);
        end
        drawimages(app);
        return;
    end

    % ---- SINGLE CLICK = START DRAG WINDOWING ----
    % Left button only if available
    if isprop(event,'Button') && event.Button ~= 1
        return;
    end

    s.isDragging = true;
    s.startPoint = app.UIFigure.CurrentPoint;
    s.startWW = app.WindowWidthSlider.Value;     % Window Width
    s.startWL = app.WindowLevelSlider.Value;   % Window Level
    app.UIFigure.UserData = s;

    figure(app.UIFigure);
        end

        % Window button up function: UIFigure
        function UIFigureWindowButtonUp(app, event)
           s = app.UIFigure.UserData;
    if isstruct(s)
        s.isDragging = false;
        app.UIFigure.UserData = s;
    end
        end

        % Window button motion function: UIFigure
        function UIFigureWindowButtonMotion(app, event)
         s = app.UIFigure.UserData;
    if ~isstruct(s) || ~isfield(s,'isDragging') || ~s.isDragging
        return;
    end

    cp = app.UIFigure.CurrentPoint;           % [x y] in figure pixels
    dx = cp(1) - s.startPoint(1);             % horizontal drag
    dy = cp(2) - s.startPoint(2);             % vertical drag

    % ---- PACS-like sensitivity ----
    k = 1/300;   % increase (e.g. 1/200) if you want faster response

    % Window Width: multiplicative (exponential) so it never goes <= 0
    newWW = s.startWW * exp(dx * k);

    % Window Level: additive, scaled by starting WW (so it feels consistent)
    newWL = s.startWL + dy * (s.startWW * k);

    % ---- Clamp to slider limits ----
    newWW = max(app.WindowWidthSlider.Limits(1),   min(app.WindowWidthSlider.Limits(2),   newWW));
    newWL = max(app.WindowLevelSlider.Limits(1), min(app.WindowLevelSlider.Limits(2), newWL));

    % ---- Update sliders (WW = ContrastSlider, WL = BrightnessSlider) ----
    app.WindowWidthSlider.Value   = newWW;
    app.WindowLevelSlider.Value = newWL;

    updateWindowOnly(app);
        end

        % Window scroll wheel function: UIFigure
        function UIFigureWindowScrollWheel(app, event)
             % No data loaded
    if isempty(app.scan_data) || isempty(app.currentscan)
        return;
    end

    % Wheel direction (+ve = scroll down)
    step = event.VerticalScrollCount;

    % Keyboard modifiers
    mods = app.UIFigure.CurrentModifier;
    shiftHeld = any(strcmp(mods, 'shift'));
    ctrlHeld  = any(strcmp(mods, 'control'));

    % Decide what we are scrolling
    if ctrlHeld
        dd = app.EvolutionFieldDropDown;   % Ctrl + wheel → field
    elseif shiftHeld
        dd = app.EvolutionTimeDropDown;    % Shift + wheel → time
    else
        dd = app.SliceDropDown;            % Wheel → slice
    end

    if isempty(dd.ItemsData)
        return;
    end

    % Current value
    v = dd.Value;
    vmin = min(dd.ItemsData);
    vmax = max(dd.ItemsData);

    % Scroll (invert sign here if direction feels wrong)
    vNew = v + step;

    % Clamp
    vNew = max(vmin, min(vmax, vNew));

    if vNew ~= v
        dd.Value = vNew;

        % Field changes must update time list
        if dd == app.EvolutionFieldDropDown
            EvolutionFieldDropDownValueChanged(app, []);
        else
            drawimages(app);
        end
    end
            
        end

        % Value changed function: CombinationDropDown
        function CombinationDropDownValueChanged(app, event)
            value = app.CombinationDropDown.Value;
            
        end

        % Button pushed function: GeometricAccButton
        function GeometricAccButtonPushed(app, event)
      % ---------------------- Your existing setup ----------------------
disp_slice = app.SliceDropDown.Value;
disp_field = app.EvolutionFieldDropDown.Value;
disp_time  = app.EvolutionTimeDropDown.Value;

targetSlice = abs(app.scan_data{app.currentscan}.magimage(:,:,disp_slice,disp_time,disp_field));

[rows, cols] = size(targetSlice);

fov = app.scan_data{app.currentscan}.fov;
pixelSize = fov/rows;

targetSliceNorm = mat2gray(targetSlice);

% ---------------------- Phantom mask (keep your logic, but add 2 safety lines) ----------------------
maskInitial = imbinarize(targetSliceNorm);
maskInitial = imfill(maskInitial,'holes');
maskInitial = bwareaopen(maskInitial, round(0.001*rows*cols));

se = strel('disk', 3);
maskFilled = imclose(maskInitial, se);

CC = bwconncomp(maskFilled);
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);

maskMROI = false(size(maskFilled));
maskMROI(CC.PixelIdxList{idx}) = true;

% Need stats early for phantom centroid fallback
phantomStats = regionprops(maskMROI,'Centroid');
cx0 = phantomStats.Centroid(1);
cy0 = phantomStats.Centroid(2);

% ---------------------- GRID CENTRE (projection method) ----------------------
I = targetSliceNorm;

maskInner = imerode(maskMROI, strel('disk', 12));   % tweak 8–20
if ~any(maskInner(:)), maskInner = maskMROI; end

Iroi = I; 
Iroi(~maskInner) = 0;

Ibh = imbothat(Iroi, strel('disk', 7));            % tweak 3–9
Ibh = mat2gray(Ibh);
Ibh = imgaussfilt(Ibh, 0.8);                       % tweak 0–1.5

kv = [1 -2 1; 1 -2 1; 1 -2 1];
kh = kv';

Ev = abs(imfilter(Ibh, kv, 'replicate'));
Eh = abs(imfilter(Ibh, kh, 'replicate'));
Ev(~maskInner) = 0;
Eh(~maskInner) = 0;

profX = smoothdata(sum(Ev,1), 'gaussian', 11);
profY = smoothdata(sum(Eh,2), 'gaussian', 11);

colMask = any(maskInner,1);
rowMask = any(maskInner,2);

[pksX, locX] = findpeaks(profX .* colMask, 'MinPeakProminence', 0.02*max(profX));
[pksY, locY] = findpeaks(profY .* rowMask, 'MinPeakProminence', 0.02*max(profY));

if isempty(locX) || isempty(locY)
    warning('Grid projections failed; falling back to phantom centroid.');
    x_center = cx0;
    y_center = cy0;
else
    bestScore = -Inf;
    bestXY = [cx0, cy0];

    distPenalty = 1; % <<<< FIXED (your "1" made it overly aggressive)

    for i = 1:numel(locX)
        for j = 1:numel(locY)
            x = locX(i);
            y = locY(j);

            xi = max(1, min(cols, round(x)));
            yi = max(1, min(rows, round(y)));

            if ~maskInner(yi, xi)
                continue
            end

            d2 = (x - cx0)^2 + (y - cy0)^2;
            score = (pksX(i) + pksY(j)) - d2/distPenalty;

            if score > bestScore
                bestScore = score;
                bestXY = [x, y];
            end
        end
    end

    x_center = bestXY(1);
    y_center = bestXY(2);
end

gridCenter = [x_center, y_center];

% ---------------------- GRID ANGLE (constrained Hough + intersection suppression) ----------------------
lineLen = 15;   % tweak 10–25 depending on resolution
horizLines = imopen(Ibh, strel('line', lineLen, 0));
vertLines  = imopen(Ibh, strel('line', lineLen, 90));
Igrid = max(horizLines, vertLines);
Igrid(~maskInner) = 0;

edges = edge(Igrid, 'canny');
edges(~maskInner) = 0;

thetaRange = [-10:0.5:10, 80:0.5:89, -90:0.5:-80];
[H,T,R] = hough(edges, 'Theta', thetaRange);

P = houghpeaks(H, 25, 'Threshold', 0.30*max(H(:)));
hlines = houghlines(edges, T, R, P, 'FillGap', 25, 'MinLength', 60);

angles = [];
weights = [];

for k = 1:numel(hlines)
    dx = hlines(k).point2(1) - hlines(k).point1(1);
    dy = hlines(k).point2(2) - hlines(k).point1(2);
    Lseg = hypot(dx, dy);
    if Lseg < 1, continue; end

    ang = atan2d(dy, dx);     % line direction
    ang = mod(ang, 90);       % fold to [0,90)

    angles(end+1)  = ang; %#ok<AGROW>
    weights(end+1) = Lseg; %#ok<AGROW>
end

if isempty(angles)
    gridAngleDeg = 0;
else
    S = sum(weights .* sind(2*angles));
    C = sum(weights .* cosd(2*angles));
    gridAngleDeg = 0.5 * atan2d(S, C);
end

% ---------------------- DIAMETER MEASUREMENT ALONG 45°-OFFSET AXES ----------------------
% Missing fluid shortens one of the grid-aligned axes when clipped to mask.
% Use axes rotated by +45° (and +135°) relative to the grid axes.
phiMeas = gridAngleDeg + 45;

L = 0.90 * min(rows, cols);     % long enough to span phantom
n = 4000;                       % sampling density

% Build two orthogonal measurement lines (unclipped)
xA = [x_center - L*cosd(phiMeas),      x_center + L*cosd(phiMeas)];
yA = [y_center - L*sind(phiMeas),      y_center + L*sind(phiMeas)];
xB = [x_center - L*cosd(phiMeas+90),   x_center + L*cosd(phiMeas+90)];
yB = [y_center - L*sind(phiMeas+90),   y_center + L*sind(phiMeas+90)];

% Clip A to maskMROI by sampling along the line
xs = linspace(xA(1), xA(2), n);
ys = linspace(yA(1), yA(2), n);
xi = max(1, min(cols, round(xs)));
yi = max(1, min(rows, round(ys)));
inside = maskMROI(sub2ind([rows, cols], yi, xi));
idx = find(inside);
if ~isempty(idx)
    i1 = idx(1); i2 = idx(end);
    xAc = [xs(i1), xs(i2)];
    yAc = [ys(i1), ys(i2)];
    diamA_mm = hypot(xAc(2)-xAc(1), yAc(2)-yAc(1)) * pixelSize;
else
    xAc = xA; yAc = yA;
    diamA_mm = NaN;
end

% Clip B to maskMROI
xs = linspace(xB(1), xB(2), n);
ys = linspace(yB(1), yB(2), n);
xi = max(1, min(cols, round(xs)));
yi = max(1, min(rows, round(ys)));
inside = maskMROI(sub2ind([rows, cols], yi, xi));
idx = find(inside);
if ~isempty(idx)
    i1 = idx(1); i2 = idx(end);
    xBc = [xs(i1), xs(i2)];
    yBc = [ys(i1), ys(i2)];
    diamB_mm = hypot(xBc(2)-xBc(1), yBc(2)-yBc(1)) * pixelSize;
else
    xBc = xB; yBc = yB;
    diamB_mm = NaN;
end

% Use the measurement that is LESS affected by the missing-fluid notch.
% Typically the "too short" one will be the minimum, so use the MAX.
calcDiameter = max(diamA_mm, diamB_mm);

% ---------------------- Plot + draw 45°-offset measurement axes ----------------------
acrDiameter = 190;

figure('Name', 'Geometric Accuracy Analysis', 'NumberTitle', 'off');
imshow(labeloverlay(targetSliceNorm, maskMROI, 'Transparency', 0.8, 'Colormap', 'spring'));
axis square;

title(sprintf('Measured Diameter (45° axes): %.2f mm  (A=%.2f, B=%.2f)  Angle=%.1f°', ...
    calcDiameter, diamA_mm, diamB_mm, phiMeas));

hold on
plot(x_center, y_center, 'r+', 'MarkerSize', 14, 'LineWidth', 2);

line(xAc, yAc, 'Color', 'cyan', 'LineWidth', 2, 'DisplayName', sprintf('Meas axis %.1f°', mod(phiMeas,180)));
line(xBc, yBc, 'Color', 'cyan', 'LineWidth', 2, 'DisplayName', sprintf('Meas axis %.1f°', mod(phiMeas+90,180)));

hold off

% ---------------------- PASS/MARGINAL/FAIL (NOW BASED ON 45° DIAMETER) ----------------------
err = abs(calcDiameter - acrDiameter);
if err <= 2
    resultStr = 'PASS';
elseif err <= 3
    resultStr = 'MARGINAL';
else
    resultStr = 'FAIL';
end

fprintf('--- Geometric Accuracy (45° axes) ---\n');
fprintf('Grid centre (px): [%.1f, %.1f]\n', x_center, y_center);
fprintf('Grid angle (deg): %.2f\n', gridAngleDeg);
fprintf('Measured axes angle (deg): %.2f\n', phiMeas);
fprintf('Diameter A (mm): %.2f\nDiameter B (mm): %.2f\n', diamA_mm, diamB_mm);
fprintf('Selected Diameter (mm): %.2f\nACR Diameter (mm): %.2f\nError (mm): %.2f\nResult: %s\n', ...
    calcDiameter, acrDiameter, err, resultStr);
        end

        % Button pushed function: UniformityButton
        function UniformityButtonPushed(app, event)
             %This can be hardcoded for the specific slice once confirmed
            disp_slice = app.SliceDropDown.Value;
            disp_field = app.EvolutionFieldDropDown.Value;
            disp_time = app.EvolutionTimeDropDown.Value;

            targetSlice = abs(app.scan_data{app.currentscan}.magimage(:,:,disp_slice,disp_time,disp_field));

            %Extracts the dimesnions fro  the scan data
            rows = app.scan_data{app.currentscan}.views;
            cols = app.scan_data{app.currentscan}.samples;

            fov = app.scan_data{app.currentscan}.fov;
            %Calculating the size of one pixel in mm
            pixelSize = fov/rows;

            %Normalising the pixel values within the slice
            targetSliceNorm = mat2gray(targetSlice);
            %Thresholding using an adaptive sensitivity
            maskInitial = imbinarize(targetSliceNorm); %Add this line of code if the mask isnt correct , "adaptive", Sensitivity=0.7);
            se = strel('disk', 3);
            %Making a smoothed mask of the signal region of the phantom
            maskFilled = imclose(maskInitial, se);

            %Identify and count pixels in connected groups of pixels in the binary mask
            CC = bwconncomp(maskFilled);
            %Applying the numel function to every cell (group of connected pixels) in the CC array
            numPixels = cellfun(@numel, CC.PixelIdxList);
            %Finding the largest region of pixels; the phantom
            [~, idx] = max(numPixels);

            %Create a binary mask with the largest region of pixels
            maskMROI = zeros(rows, cols, 'logical');
            maskMROI(CC.PixelIdxList{idx}) = true;

            %Shrink mask by 20mm from the phantom edge
            shrinkPixels = round(20 / pixelSize);
            maskMROIRefined = imerode(maskMROI, strel('disk', shrinkPixels));
           
            %Create a circular averaging filter to simulate a 1cm^2 ROI (ACR guidance)
            %Get the highest region rather than single pixel to avoid noise outliers
            roiSize = round(sqrt(100/pi) / pixelSize);
            %Making a smoothing filter, disk shaped mask with the specified radius
            h = fspecial('disk', roiSize);
            %imfilter moves the mask over every pixel and replaces the value with the average of its neighbours in the given radius
            %'replicate' assumes the edge pixels of the image are the same as the overhang
            smoothedSlice = imfilter(targetSlice, h, 'replicate');

            %Find high and low mean signal values within the smoothed signal region
            phantomValues = smoothedSlice(maskMROIRefined);
            high = max(phantomValues);
            low = min(phantomValues);

            %Locate the centers of these regions for plotting, by finding pixel equal to the high or low that falls in the refined mask region
            [maxRow, maxCol] = find(smoothedSlice == high & maskMROIRefined, 1);
            [minRow, minCol] = find(smoothedSlice == low & maskMROIRefined, 1);

            %Percent Integral Uniformity (ACR guidance)
            piu = 100 * (1 - ((high - low)/(high + low)));

            figure('Name', 'Image Uniformity Analysis', 'NumberTitle', 'off');
            %Displaying the original slice with the refined mask overlay
            imshow(labeloverlay(targetSliceNorm, maskMROI, 'Transparency', 0.8, 'Colormap', 'summer'));
            title(sprintf('PIU: %.2f%% (High: %.2f, Low: %.2f)', piu, high, low));
            axis square;
            hold on;
            %creating a circular element which had the radius required, to
            %display onto the figure
            viscircles([maxCol, maxRow], roiSize, 'EdgeColor', 'r', 'LineWidth', 2);
            viscircles([minCol, minRow], roiSize, 'EdgeColor', 'b', 'LineWidth', 2);
            %Non numerical 'placeholder' for the legend, doesn't work with
            %viscircles
            hMax = plot(NaN, NaN, 'or', 'LineWidth', 2);
            hMin = plot(NaN, NaN, 'ob', 'LineWidth', 2); 
            hMask = plot(NaN, NaN, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 1 0.8]);
            legend([hMask, hMax, hMin], {'Refined ROI', 'Max Signal ROI', 'Min Signal ROI'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
            hold off;

            %For images acquired using the large phantom on MRI systems with field strengths less than 3 Tesla, PIU should be greater than or equal to 87.5% and will fail if less than 85%.
            fprintf('--- Image Uniformity ---\n');
            fprintf('PIU: %.2f%%\nMax ROI Mean: %.2f\nMin ROI Mean: %.2f\n', piu, high, low);
            
            if piu >= 87.5
                fprintf('Result: PASS\n');
            elseif piu >= 85
                fprintf('Result: MARGINAL (Pass >= 85, but below 87.5 preferred)\n');
            else
                fprintf('Result: FAIL (below 85)\n');
            end
        end

        % Button pushed function: THKButton
        function THKButtonPushed(app, event)
             %This can be hardcoded for the specific slice once confirmed
            disp_slice = app.SliceDropDown.Value;
            disp_field = app.EvolutionFieldDropDown.Value;
            disp_time = app.EvolutionTimeDropDown.Value;

            targetSlice = abs(app.scan_data{app.currentscan}.magimage(:,:,disp_slice,disp_time,disp_field));

            %Extracts the dimesnions fro  the scan data
            rows = app.scan_data{app.currentscan}.views;
            cols = app.scan_data{app.currentscan}.samples;

            fov = app.scan_data{app.currentscan}.fov;
            %Calculating the size of one pixel in mm
            pixelSize = fov/rows;

            %Normalising the pixel values within the slice
            targetSliceNorm = mat2gray(targetSlice);

            figure('Name', 'Slice Thickness - Manual')
            imagesc(targetSliceNorm)
            title('Slice Thickness - Manual',  'FontSize', 20)
            axis square
            axis off
            colormap gray
            annotation('textbox', [0.28, 0.0, 0.1, 0.1], 'String', "Draw 2 lines, measuring the lengths of the signal bars")

            length1 = drawline('Color', 'b');
            length2 = drawline('Color', 'g');

            %The position of the roi line is given by a 2x2 matrix where each row
            %is one of the end points of the line
            pos1 = length1.Position;
            pos2 = length2.Position;

            %Use the geometric distance formula to calculate the length in pixels
            pixelLength1 = sqrt(diff(pos1(:,1))^2 + diff(pos1(:,2))^2);
            pixelLength2 = sqrt(diff(pos2(:,1))^2 + diff(pos2(:,2))^2);

            top = pixelSize*pixelLength1;
            bottom = pixelSize*pixelLength2;

            %Calculate the slice thickness
            sliceThickness = 0.2 * ( (top * bottom) / (top + bottom) );
            trueThickness = 5.0; %mm

            %Print the result
            fprintf('Slice Thickness: %.2f mm\n', sliceThickness);

            %ACR guidance states that all measured lengths should be 5.0 mm +/- 0.7 mm.
            %It will fail if it is +/- 1.0 mm
            if abs(sliceThickness - trueThickness) <= 0.7
                fprintf('Result: PASS\n');
            elseif abs(sliceThickness - trueThickness) >= 0.7 && abs(sliceThickness - trueThickness) <= 1.0
                fprintf('Result: Over 0.7 mm difference, but under 1.0 mm failure point\n');
            else
                fprintf('Result: FAIL (Outwith 1 mm tolerance)\n');
            end
        end

        % Button pushed function: AllButton
        function AllButtonPushed(app, event)
            %DEFINITIONS
            %Current GUI image
            disp_slice = app.SliceDropDown.Value;
            disp_field = app.EvolutionFieldDropDown.Value;
            disp_time = app.EvolutionTimeDropDown.Value;

            %Hardcoded slices
            sliceSNR = 41;
            sliceGEO = 41;
            sliceUNI = 41;
            sliceTHI = 70;

            targetSliceSNR = abs(app.scan_data{app.currentscan}.magimage(:,:,sliceSNR,disp_time,disp_field));
            targetSliceGEO = abs(app.scan_data{app.currentscan}.magimage(:,:,sliceGEO,disp_time,disp_field));
            targetSliceUNI = abs(app.scan_data{app.currentscan}.magimage(:,:,sliceUNI,disp_time,disp_field));
            targetSliceTHI = abs(app.scan_data{app.currentscan}.magimage(:,:,sliceTHI,disp_time,disp_field));

            targetSliceSNRNorm = mat2gray(targetSliceSNR);
            targetSliceGEONorm = mat2gray(targetSliceGEO);
            targetSliceUNINorm = mat2gray(targetSliceUNI);
            targetSliceTHINorm = mat2gray(targetSliceTHI);

            %Dimensions of the scan data
            rows = app.scan_data{app.currentscan}.views;
            cols = app.scan_data{app.currentscan}.samples;

            fov = app.scan_data{app.currentscan}.fov;

            pixelSize = fov/rows;

            acrDiameter = 190;

%PRE FIGURE CALCULATIONS
            maskInitialSNR = imbinarize(targetSliceSNRNorm); %Add this line of code if the mask isnt correct , "adaptive", Sensitivity=0.7);
            seSNR = strel('disk', 4);
            maskFilledSNR = imclose(maskInitialSNR,seSNR);
            
            %Identify and count pixels in connected groups of pixels in the binary mask
            CCSNR = bwconncomp(maskFilledSNR);
            %Applying the numel function to every cell (group of connected pixels) in
            %the CC array
            numPixelsSNR = cellfun(@numel, CCSNR.PixelIdxList);
            %Finding the largest region of pixels; the phantom
            [~, idx] = max(numPixelsSNR);

            %Create a binary mask with the largest region of pixels
            maskMROISNR = zeros(rows, cols, 'logical');
            maskMROISNR(CCSNR.PixelIdxList{idx}) = true;
    
            %Using imerode to shrink the mask to ensure none of the edge pixels are included in the signal average
            seErodeSNR = strel('disk', 6);
            maskMROIRefinedSNR = imerode(maskMROISNR, seErodeSNR);
            
            %Mean signal intensity
            S = mean(targetSliceSNR(maskMROIRefinedSNR));
    
            %Number and size of noise ROI's
            noiseRegions = 4;
            noiseStd = zeros(noiseRegions, 1);
            noiseBoxSize = round(rows * 0.05);
            buffer = 8;
    
            %Creating a table of starting coordinates to place the noise regions, with an 8 pixel buffer from the edges
            cornerCoords = [
                buffer, buffer;
                cols - noiseBoxSize - buffer, buffer;
                cols - noiseBoxSize - buffer, rows - noiseBoxSize - buffer;
                buffer, rows - noiseBoxSize - buffer
            ];

            %Applying a more sensitive threshold
            maskInitialGEO = imbinarize(targetSliceGEONorm); %Add this line of code if the mask isnt correct , "adaptive", Sensitivity=0.7);
            seGEO = strel('disk', 3);
            maskFilledGEO = imclose(maskInitialGEO, seGEO);

            %Identify and count pixels in connected groups of pixels in the binary mask
            CCGEO = bwconncomp(maskFilledGEO);
            %Applying the numel function to every cell (group of connected pixels) in the CC array
            numPixelsGEO = cellfun(@numel, CCGEO.PixelIdxList);
            %Finding the largest region of pixels; the phantom
            [~, idx] = max(numPixelsGEO);

            %Create a binary mask with the largest region of pixels
            maskMROIGEO = zeros(rows, cols, 'logical');
            maskMROIGEO(CCGEO.PixelIdxList{idx}) = true;

            stats = regionprops(maskMROIGEO, 'EquivDiameter', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Centroid');
            %The diameter if the region was a perfect circle
            calcDiameter = pixelSize * stats.EquivDiameter;
            %Max and min diameters to help identify any distortion
            maxDiam = stats.MajorAxisLength * pixelSize;
            minDiam = stats.MinorAxisLength * pixelSize;
            distortion = stats.MajorAxisLength / stats.MinorAxisLength;

            %Thresholding using an adaptive sensitivity
            maskInitialUNI = imbinarize(targetSliceUNINorm); %Add this line of code if the mask isnt correct , "adaptive", Sensitivity=0.7);
            seUNI = strel('disk', 3);
            %Making a smoothed mask of the signal region of the phantom
            maskFilledUNI = imclose(maskInitialUNI, seUNI);

            %Identify and count pixels in connected groups of pixels in the binary mask
            CCUNI = bwconncomp(maskFilledUNI);
            %Applying the numel function to every cell (group of connected pixels) in the CC array
            numPixelsUNI = cellfun(@numel, CCUNI.PixelIdxList);
            %Finding the largest region of pixels; the phantom
            [~, idx] = max(numPixelsUNI);

            %Create a binary mask with the largest region of pixels
            maskMROIUNI = zeros(rows, cols, 'logical');
            maskMROIUNI(CCUNI.PixelIdxList{idx}) = true;

            %Shrink mask by 20mm from the phantom edge
            shrinkPixelsUNI = round(20 / pixelSize);
            maskMROIRefinedUNI = imerode(maskMROIUNI, strel('disk', shrinkPixelsUNI));
           
            %Create a circular averaging filter to simulate a 1cm^2 ROI (ACR guidance)
            %Get the highest region rather than single pixel to avoid noise outliers
            roiSize = round(sqrt(100/pi) / pixelSize);
            %Making a smoothing filter, disk shaped mask with the specified radius
            h = fspecial('disk', roiSize);
            %imfilter moves the mask over every pixel and replaces the value with the average of its neighbours in the given radius
            %'replicate' assumes the edge pixels of the image are the same as the overhang
            smoothedSlice = imfilter(targetSliceUNI, h, 'replicate');

            %Find high and low mean signal values within the smoothed signal region
            phantomValues = smoothedSlice(maskMROIRefinedUNI);
            high = max(phantomValues);
            low = min(phantomValues);

            %Locate the centers of these regions for plotting, by finding pixel equal to the high or low that falls in the refined mask region
            [maxRow, maxCol] = find(smoothedSlice == high & maskMROIRefinedUNI, 1);
            [minRow, minCol] = find(smoothedSlice == low & maskMROIRefinedUNI, 1);

%FIGURES
            name = sprintf('QA Test (%s)', datetime("now"));
            figure("Name", name)
            t = tiledlayout(2,2, "TileSpacing", "tight");
            axis square
            axis off
            colormap gray
            title(t,'Daily QA Test')

            %SNR
            nexttile
            imshow(labeloverlay(targetSliceSNRNorm, maskMROIRefinedSNR, 'Transparency', 0.7, 'Colormap', 'hot'));
            for i = 1:noiseRegions
                xStart = cornerCoords(i, 1);
                yStart = cornerCoords(i, 2);

                maskNoise = false(rows, cols);
                maskNoise(yStart : yStart + noiseBoxSize, xStart : xStart + noiseBoxSize) = true;
                noiseStd(i) = std(targetSliceSNR(maskNoise));

                rectangle('Position', [xStart, yStart, noiseBoxSize, noiseBoxSize], 'EdgeColor', 'g', 'LineWidth', 2);
            end
            title('SNR')

            %GEOMETRIC ACCURACY
            nexttile
            imshow(labeloverlay(targetSliceGEONorm, maskMROIGEO, 'Transparency', 0.8, 'Colormap', 'spring'));
            hold on;
            phi = -stats.Orientation;
            x_center = stats.Centroid(1);
            y_center = stats.Centroid(2);
            
            x_line = [x_center - (stats.MajorAxisLength/2) * cosd(phi), x_center + (stats.MajorAxisLength/2) * cosd(phi)];
            y_line = [y_center - (stats.MajorAxisLength/2) * sind(phi), y_center + (stats.MajorAxisLength/2) * sind(phi)];
            
            plot(stats.Centroid(1), stats.Centroid(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
            line(x_line, y_line, 'Color', 'cyan', 'LineWidth', 2, 'DisplayName', 'Major Axis');
            hold off;
            title('Geometric Accuracy')

            %IMAGE UNIFORMITY
            nexttile
            imshow(labeloverlay(targetSliceUNINorm, maskMROIUNI, 'Transparency', 0.8, 'Colormap', 'summer'));
            hold on;
            %creating a circular element which had the radius required, to display onto the figure
            viscircles([maxCol, maxRow], roiSize, 'EdgeColor', 'r', 'LineWidth', 2);
            viscircles([minCol, minRow], roiSize, 'EdgeColor', 'b', 'LineWidth', 2);
            %Non numerical 'placeholder' for the legend, doesn't work with viscircles
            hMax = plot(NaN, NaN, 'or', 'LineWidth', 2);
            hMin = plot(NaN, NaN, 'ob', 'LineWidth', 2); 
            hMask = plot(NaN, NaN, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 1 0.8]);
            legend([hMask, hMax, hMin], {'Refined ROI', 'Max Signal ROI', 'Min Signal ROI'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
            hold off;
            title('Image Uniformity')

            %SLICE THICKNESS ACCURACY
            nexttile
            imshow(targetSliceTHINorm)
            title('Slice Thickness Accuracy')
            length1 = drawline('Color', 'b');
            length2 = drawline('Color', 'g');
            

%POST FIGURE CALCULATIONS
            %Calculation of the mean noise which is corrected by dividing by 0.66 which accounts for the Rayleigh distribution of the noise in the magnitude image.
            noiseAvg = mean(noiseStd);
            noiseCorrected = noiseAvg / 0.66; 
            %SNR calculation
            snr = S / noiseCorrected;

            %Calculate Percent Integral Uniformity
            piu = 100 * (1 - ((high - low)/(high + low)));

            %The position of the roi line is given by a 2x2 matrix where each row
            %is one of the end points of the line
            pos1 = length1.Position;
            pos2 = length2.Position;

            %Use the geometric distance formula to calculate the length in pixels
            pixelLength1 = sqrt(diff(pos1(:,1))^2 + diff(pos1(:,2))^2);
            pixelLength2 = sqrt(diff(pos2(:,1))^2 + diff(pos2(:,2))^2);

            top = pixelSize*pixelLength1;
            bottom = pixelSize*pixelLength2;

            %Calculate the Slice Thickness
            sliceThickness = 0.2 * ( (top * bottom) / (top + bottom) );
            trueThickness = 5.0; %mm

%OUTPUTS
            fprintf('--- Signal to Noise Ratio ---\n');
            fprintf('SNR: %.3f | Mean Signal: %.4f\n | Mean Noise: %.4f\n', snr, S, noiseAvg);
            resSNR = snr;

            %ACR guidance states that all measured lengths should be within 2 mm of their true values and will fail if the length differs by more than 3 mm
            fprintf('--- Geometric Accuracy ---\n');
            fprintf('Equivalent Diameter: %.2f mm\nMajor Axis: %.2f mm\nMinor Axis: %.2f mm\nDistortion: %.3f\n', calcDiameter, maxDiam, minDiam, distortion);
            if abs(calcDiameter - acrDiameter) <= 2
                resGEO = "pass";
                fprintf('Result: PASS\n');
            elseif abs(calcDiameter - acrDiameter) <= 3
                resGEO = "marginal";
                fprintf('Result: MARGINAL\n');
            else
                resGEO = "fail";
                fprintf('Result: FAIL (over 3 mm)\n');
            end

            %For images acquired using the large phantom on MRI systems with field strengths less than 3 Tesla, PIU should be greater than or equal to 87.5% and will fail if less than 85%.
            fprintf('--- Image Uniformity ---\n');
            fprintf('PIU: %.2f%%\nMax ROI Mean: %.2f\nMin ROI Mean: %.2f\n', piu, high, low);
            
            if piu >= 87.5
                resUNI = "pass";
                fprintf('Result: PASS\n');
            elseif piu >= 85
                resUNI = "marginal";
                fprintf('Result: MARGINAL (Pass >= 85, but below 87.5 preferred)\n');
            else
                resUNI = "fail";
                fprintf('Result: FAIL (below 85)\n');
            end

            %ACR guidance states that all measured thicknesses should be 5.0 mm +/- 0.7 mm.
            %It will fail if it is +/- 1.0 mm
            fprintf('--- Slice Thickness Accuracy ---\n');
            fprintf('Slice Thickness: %.2f mm\n', sliceThickness);
            if abs(sliceThickness - trueThickness) <= 0.7
                resTHI = "pass";
                fprintf('Result: PASS\n');
            elseif abs(sliceThickness - trueThickness) >= 0.7 && abs(sliceThickness - trueThickness) <= 1.0
                resTHI = "marginal";
                fprintf('Result: Over 0.7 mm difference, but under 1.0 mm failure point\n');
            else
                resTHI = "fail";
                fprintf('Result: FAIL (Outwith 1 mm tolerance)\n');
            end

            results = [resSNR resGEO resUNI resTHI];
          %  writematrix(results,"C:\Users\power\MRI Project\results.xlsx", 'WriteMode', 'append');
        end

        % Button pushed function: FilterSettingsButton
        function FilterSettingsButtonPushed(app, event)
     if isempty(app.scan_data) || isempty(app.currentscan)
        return
    end

    openFilterSubGui(app);
        end

        % Value changed function: FFTDropDown_3D
        function FFTDropDown_3DValueChanged(app, event)
            updateimages(app);
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1274 746];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, @UIFigureWindowButtonDown, true);
            app.UIFigure.WindowButtonUpFcn = createCallbackFcn(app, @UIFigureWindowButtonUp, true);
            app.UIFigure.WindowButtonMotionFcn = createCallbackFcn(app, @UIFigureWindowButtonMotion, true);
            app.UIFigure.WindowScrollWheelFcn = createCallbackFcn(app, @UIFigureWindowScrollWheel, true);
            app.UIFigure.WindowKeyPressFcn = createCallbackFcn(app, @UIFigureWindowKeyPress, true);
            app.UIFigure.KeyPressFcn = [];

            % Create axis_kspace
            app.axis_kspace = uiaxes(app.UIFigure);
            app.axis_kspace.PlotBoxAspectRatio = [1 1.0020366598778 1];
            app.axis_kspace.XColor = 'none';
            app.axis_kspace.XTick = [];
            app.axis_kspace.YColor = 'none';
            app.axis_kspace.YTick = [];
            app.axis_kspace.Color = 'none';
            app.axis_kspace.ClippingStyle = 'rectangle';
            app.axis_kspace.Clipping = 'off';
            app.axis_kspace.Tag = 'axis_kspace';
            app.axis_kspace.Position = [118 347 379 380];

            % Create axis_image
            app.axis_image = uiaxes(app.UIFigure);
            app.axis_image.AmbientLightColor = 'none';
            app.axis_image.PlotBoxAspectRatio = [1 1.0020366598778 1];
            app.axis_image.XColor = 'none';
            app.axis_image.XTick = [];
            app.axis_image.YColor = 'none';
            app.axis_image.YTick = [];
            app.axis_image.Color = 'none';
            app.axis_image.ClippingStyle = 'rectangle';
            app.axis_image.Clipping = 'off';
            app.axis_image.ButtonDownFcn = createCallbackFcn(app, @axis_imageButtonDown, true);
            app.axis_image.Tag = 'axis_image';
            app.axis_image.Position = [533 347 379 380];

            % Create axis_coordinates
            app.axis_coordinates = uiaxes(app.UIFigure);
            app.axis_coordinates.AmbientLightColor = 'none';
            app.axis_coordinates.PlotBoxAspectRatio = [1 1.0020366598778 1];
            app.axis_coordinates.XColor = 'none';
            app.axis_coordinates.XTick = [];
            app.axis_coordinates.YColor = 'none';
            app.axis_coordinates.YTick = [];
            app.axis_coordinates.Color = 'none';
            app.axis_coordinates.ClippingStyle = 'rectangle';
            app.axis_coordinates.Tag = 'axis_image';
            app.axis_coordinates.Position = [921 655 67 72];

            % Create panel_recon
            app.panel_recon = uipanel(app.UIFigure);
            app.panel_recon.Title = 'Reconstruction Options';
            app.panel_recon.Tag = 'panel_recon';
            app.panel_recon.Scrollable = 'on';
            app.panel_recon.Position = [25 10 174 327];

            % Create Sizeof2DFFTLabel
            app.Sizeof2DFFTLabel = uilabel(app.panel_recon);
            app.Sizeof2DFFTLabel.Enable = 'off';
            app.Sizeof2DFFTLabel.Position = [13 265 46 22];
            app.Sizeof2DFFTLabel.Text = '2D FFT';

            % Create FFTDropDown
            app.FFTDropDown = uidropdown(app.panel_recon);
            app.FFTDropDown.Items = {'64', '128', '256', '512'};
            app.FFTDropDown.ValueChangedFcn = createCallbackFcn(app, @FFTDropDownValueChanged, true);
            app.FFTDropDown.Enable = 'off';
            app.FFTDropDown.Tooltip = {'Size of the FFT. Larger sizes of transform are accomplished by symmetrically zero-padding kspace'};
            app.FFTDropDown.Position = [95 265 67 22];
            app.FFTDropDown.Value = '64';

            % Create WindowDropDownLabel
            app.WindowDropDownLabel = uilabel(app.panel_recon);
            app.WindowDropDownLabel.Enable = 'off';
            app.WindowDropDownLabel.Position = [13 197 48 22];
            app.WindowDropDownLabel.Text = 'Window';

            % Create WindowDropDown
            app.WindowDropDown = uidropdown(app.panel_recon);
            app.WindowDropDown.Items = {'None', 'Hann', 'Blackman', 'Hamming', 'Kaiser-Bessel', 'Tukey'};
            app.WindowDropDown.ValueChangedFcn = createCallbackFcn(app, @WindowDropDownValueChanged, true);
            app.WindowDropDown.Enable = 'off';
            app.WindowDropDown.Tooltip = {'kspace windowing function'};
            app.WindowDropDown.Position = [95 197 67 22];
            app.WindowDropDown.Value = 'None';

            % Create WindowSizeEditFieldLabel
            app.WindowSizeEditFieldLabel = uilabel(app.panel_recon);
            app.WindowSizeEditFieldLabel.Enable = 'off';
            app.WindowSizeEditFieldLabel.Position = [13 163 75 22];
            app.WindowSizeEditFieldLabel.Text = 'Window Size';

            % Create WindowSizeEditField
            app.WindowSizeEditField = uieditfield(app.panel_recon, 'numeric');
            app.WindowSizeEditField.Limits = [0 10];
            app.WindowSizeEditField.ValueChangedFcn = createCallbackFcn(app, @WindowSizeEditFieldValueChanged, true);
            app.WindowSizeEditField.HorizontalAlignment = 'center';
            app.WindowSizeEditField.Enable = 'off';
            app.WindowSizeEditField.Tooltip = {'Size of the windowing function as a fraction of the base image size'};
            app.WindowSizeEditField.Position = [95 163 67 22];
            app.WindowSizeEditField.Value = 1;

            % Create RebuildImagesButton
            app.RebuildImagesButton = uibutton(app.panel_recon, 'push');
            app.RebuildImagesButton.ButtonPushedFcn = createCallbackFcn(app, @RebuildImagesButtonPushed, true);
            app.RebuildImagesButton.Enable = 'off';
            app.RebuildImagesButton.Tooltip = {'Rebuilds images from phase-corrected k-space. Useful when testing new reconstruction methods'};
            app.RebuildImagesButton.Position = [11 9 145 22];
            app.RebuildImagesButton.Text = 'Rebuild Images';

            % Create CombinationDropDownLabel
            app.CombinationDropDownLabel = uilabel(app.panel_recon);
            app.CombinationDropDownLabel.Enable = 'off';
            app.CombinationDropDownLabel.Position = [13 129 72 22];
            app.CombinationDropDownLabel.Text = 'Combination';

            % Create CombinationDropDown
            app.CombinationDropDown = uidropdown(app.panel_recon);
            app.CombinationDropDown.Items = {'RSSQ', 'Adaptive', 'ACS/Walsh'};
            app.CombinationDropDown.ValueChangedFcn = createCallbackFcn(app, @CombinationDropDownValueChanged, true);
            app.CombinationDropDown.Enable = 'off';
            app.CombinationDropDown.Tooltip = {'Algorithm for combining multiple receivers'};
            app.CombinationDropDown.Position = [95 129 67 22];
            app.CombinationDropDown.Value = 'Adaptive';

            % Create BiasFieldDropDownLabel
            app.BiasFieldDropDownLabel = uilabel(app.panel_recon);
            app.BiasFieldDropDownLabel.Enable = 'off';
            app.BiasFieldDropDownLabel.Position = [13 95 61 22];
            app.BiasFieldDropDownLabel.Text = 'Bias Field';

            % Create BiasFieldDropDown
            app.BiasFieldDropDown = uidropdown(app.panel_recon);
            app.BiasFieldDropDown.Items = {'Off', 'ACS'};
            app.BiasFieldDropDown.Enable = 'off';
            app.BiasFieldDropDown.Tooltip = {'Attempt to correct uniformity inhomogeneity due to coil sensitivities'};
            app.BiasFieldDropDown.Position = [95 94 67 22];
            app.BiasFieldDropDown.Value = 'Off';

            % Create FilterSettingsButton
            app.FilterSettingsButton = uibutton(app.panel_recon, 'push');
            app.FilterSettingsButton.ButtonPushedFcn = createCallbackFcn(app, @FilterSettingsButtonPushed, true);
            app.FilterSettingsButton.Enable = 'off';
            app.FilterSettingsButton.Tooltip = {'Rebuilds images from phase-corrected k-space. Useful when testing new reconstruction methods'};
            app.FilterSettingsButton.Position = [11 41 145 23];
            app.FilterSettingsButton.Text = 'Denoising';

            % Create Sizeof3DFFTLabel
            app.Sizeof3DFFTLabel = uilabel(app.panel_recon);
            app.Sizeof3DFFTLabel.Enable = 'off';
            app.Sizeof3DFFTLabel.Position = [13 231 46 22];
            app.Sizeof3DFFTLabel.Text = '3D FFT';

            % Create FFTDropDown_3D
            app.FFTDropDown_3D = uidropdown(app.panel_recon);
            app.FFTDropDown_3D.Items = {'64', '128', '256', '512'};
            app.FFTDropDown_3D.ValueChangedFcn = createCallbackFcn(app, @FFTDropDown_3DValueChanged, true);
            app.FFTDropDown_3D.Enable = 'off';
            app.FFTDropDown_3D.Tooltip = {'Size of the FFT. Larger sizes of transform are accomplished by symmetrically zero-padding kspace'};
            app.FFTDropDown_3D.Position = [95 231 67 22];
            app.FFTDropDown_3D.Value = '64';

            % Create CurrentImagePanel
            app.CurrentImagePanel = uipanel(app.UIFigure);
            app.CurrentImagePanel.Title = 'Current Image';
            app.CurrentImagePanel.Scrollable = 'on';
            app.CurrentImagePanel.Position = [217 210 418 127];

            % Create EvolutionFieldDropDownLabel
            app.EvolutionFieldDropDownLabel = uilabel(app.CurrentImagePanel);
            app.EvolutionFieldDropDownLabel.HorizontalAlignment = 'right';
            app.EvolutionFieldDropDownLabel.Enable = 'off';
            app.EvolutionFieldDropDownLabel.Position = [9 75 84 22];
            app.EvolutionFieldDropDownLabel.Text = 'Evolution Field';

            % Create EvolutionFieldDropDown
            app.EvolutionFieldDropDown = uidropdown(app.CurrentImagePanel);
            app.EvolutionFieldDropDown.Items = {'200'};
            app.EvolutionFieldDropDown.ValueChangedFcn = createCallbackFcn(app, @EvolutionFieldDropDownValueChanged, true);
            app.EvolutionFieldDropDown.Enable = 'off';
            app.EvolutionFieldDropDown.Tooltip = {'Evolution field in mT'};
            app.EvolutionFieldDropDown.Position = [108 75 79 22];
            app.EvolutionFieldDropDown.Value = '200';

            % Create EvolutionTimeDropDownLabel
            app.EvolutionTimeDropDownLabel = uilabel(app.CurrentImagePanel);
            app.EvolutionTimeDropDownLabel.HorizontalAlignment = 'right';
            app.EvolutionTimeDropDownLabel.Enable = 'off';
            app.EvolutionTimeDropDownLabel.Position = [10 45 84 22];
            app.EvolutionTimeDropDownLabel.Text = 'Evolution Time';

            % Create EvolutionTimeDropDown
            app.EvolutionTimeDropDown = uidropdown(app.CurrentImagePanel);
            app.EvolutionTimeDropDown.Items = {'0'};
            app.EvolutionTimeDropDown.ValueChangedFcn = createCallbackFcn(app, @EvolutionTimeDropDownValueChanged, true);
            app.EvolutionTimeDropDown.Enable = 'off';
            app.EvolutionTimeDropDown.Tooltip = {'Evolution time in ms'};
            app.EvolutionTimeDropDown.Position = [109 45 78 22];
            app.EvolutionTimeDropDown.Value = '0';

            % Create SliceDropDownLabel
            app.SliceDropDownLabel = uilabel(app.CurrentImagePanel);
            app.SliceDropDownLabel.HorizontalAlignment = 'right';
            app.SliceDropDownLabel.Enable = 'off';
            app.SliceDropDownLabel.Position = [61 15 32 22];
            app.SliceDropDownLabel.Text = 'Slice';

            % Create SliceDropDown
            app.SliceDropDown = uidropdown(app.CurrentImagePanel);
            app.SliceDropDown.Items = {'1'};
            app.SliceDropDown.ValueChangedFcn = createCallbackFcn(app, @SliceDropDownValueChanged, true);
            app.SliceDropDown.Enable = 'off';
            app.SliceDropDown.Tooltip = {'Slice index'};
            app.SliceDropDown.Position = [108 15 79 22];
            app.SliceDropDown.Value = '1';

            % Create ImageTypeDropDownLabel
            app.ImageTypeDropDownLabel = uilabel(app.CurrentImagePanel);
            app.ImageTypeDropDownLabel.HorizontalAlignment = 'right';
            app.ImageTypeDropDownLabel.Enable = 'off';
            app.ImageTypeDropDownLabel.Position = [210 75 68 22];
            app.ImageTypeDropDownLabel.Text = 'Image Type';

            % Create ImageTypeDropDown
            app.ImageTypeDropDown = uidropdown(app.CurrentImagePanel);
            app.ImageTypeDropDown.Items = {'Magnitude', 'Phase', 'Relaxometry Map'};
            app.ImageTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @ImageTypeDropDownValueChanged, true);
            app.ImageTypeDropDown.Enable = 'off';
            app.ImageTypeDropDown.Tooltip = {'Displayed image can be the magnitude image, a phase map or a T1 or R1 parameter map'};
            app.ImageTypeDropDown.Position = [293 75 100 22];
            app.ImageTypeDropDown.Value = 'Magnitude';

            % Create RegisterImagesButton
            app.RegisterImagesButton = uibutton(app.CurrentImagePanel, 'push');
            app.RegisterImagesButton.ButtonPushedFcn = createCallbackFcn(app, @RegisterImagesButtonPushed, true);
            app.RegisterImagesButton.Enable = 'off';
            app.RegisterImagesButton.Position = [293 47 103 22];
            app.RegisterImagesButton.Text = 'Register Images';

            % Create DisplayOptionsPanel
            app.DisplayOptionsPanel = uipanel(app.UIFigure);
            app.DisplayOptionsPanel.Title = 'Display Options';
            app.DisplayOptionsPanel.Scrollable = 'on';
            app.DisplayOptionsPanel.Position = [217 10 589 109];

            % Create WindowLevelSliderLabel
            app.WindowLevelSliderLabel = uilabel(app.DisplayOptionsPanel);
            app.WindowLevelSliderLabel.HorizontalAlignment = 'right';
            app.WindowLevelSliderLabel.Enable = 'off';
            app.WindowLevelSliderLabel.Position = [10 63 80 22];
            app.WindowLevelSliderLabel.Text = 'Window Level';

            % Create WindowLevelSlider
            app.WindowLevelSlider = uislider(app.DisplayOptionsPanel);
            app.WindowLevelSlider.Limits = [-1 1];
            app.WindowLevelSlider.MajorTicks = [];
            app.WindowLevelSlider.MajorTickLabels = {''};
            app.WindowLevelSlider.ValueChangedFcn = createCallbackFcn(app, @WindowLevelSliderValueChanged, true);
            app.WindowLevelSlider.MinorTicks = [-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
            app.WindowLevelSlider.Enable = 'off';
            app.WindowLevelSlider.Position = [97 72 189 3];

            % Create WindowWidthSliderLabel
            app.WindowWidthSliderLabel = uilabel(app.DisplayOptionsPanel);
            app.WindowWidthSliderLabel.HorizontalAlignment = 'right';
            app.WindowWidthSliderLabel.Enable = 'off';
            app.WindowWidthSliderLabel.Position = [9 35 82 22];
            app.WindowWidthSliderLabel.Text = 'Window Width';

            % Create WindowWidthSlider
            app.WindowWidthSlider = uislider(app.DisplayOptionsPanel);
            app.WindowWidthSlider.Limits = [0 1];
            app.WindowWidthSlider.MajorTicks = [];
            app.WindowWidthSlider.MajorTickLabels = {''};
            app.WindowWidthSlider.ValueChangedFcn = createCallbackFcn(app, @WindowWidthSliderValueChanged, true);
            app.WindowWidthSlider.MinorTicks = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1];
            app.WindowWidthSlider.Enable = 'off';
            app.WindowWidthSlider.Position = [109 44 177 3];
            app.WindowWidthSlider.Value = 1;

            % Create RotateImageButton
            app.RotateImageButton = uibutton(app.DisplayOptionsPanel, 'push');
            app.RotateImageButton.ButtonPushedFcn = createCallbackFcn(app, @RotateImageButtonPushed, true);
            app.RotateImageButton.Enable = 'off';
            app.RotateImageButton.Position = [318 63 100 22];
            app.RotateImageButton.Text = 'Rotate Image';

            % Create FlipHorizontallyButton
            app.FlipHorizontallyButton = uibutton(app.DisplayOptionsPanel, 'push');
            app.FlipHorizontallyButton.ButtonPushedFcn = createCallbackFcn(app, @FlipHorizontallyButtonPushed, true);
            app.FlipHorizontallyButton.Enable = 'off';
            app.FlipHorizontallyButton.Position = [318 35 101 22];
            app.FlipHorizontallyButton.Text = 'Flip Horizontally';

            % Create FlipVerticallyButton
            app.FlipVerticallyButton = uibutton(app.DisplayOptionsPanel, 'push');
            app.FlipVerticallyButton.ButtonPushedFcn = createCallbackFcn(app, @FlipVerticallyButtonPushed, true);
            app.FlipVerticallyButton.Enable = 'off';
            app.FlipVerticallyButton.Position = [318 8 100 22];
            app.FlipVerticallyButton.Text = 'Flip Vertically';

            % Create IndividualNormalisationButton
            app.IndividualNormalisationButton = uibutton(app.DisplayOptionsPanel, 'state');
            app.IndividualNormalisationButton.ValueChangedFcn = createCallbackFcn(app, @IndividualNormalisationButtonValueChanged, true);
            app.IndividualNormalisationButton.Tooltip = {'Displayed images can either be individually normalised (using imagesc) or else scaled to the maximum intensity of the first image in the data set (usually an image acquired at 200 mt)'};
            app.IndividualNormalisationButton.Enable = 'off';
            app.IndividualNormalisationButton.Text = 'Individual Normalisation';
            app.IndividualNormalisationButton.Position = [434 63 143 22];

            % Create ScanListPanel
            app.ScanListPanel = uipanel(app.UIFigure);
            app.ScanListPanel.Title = 'Scan List';
            app.ScanListPanel.Scrollable = 'on';
            app.ScanListPanel.Position = [956 336 299 306];

            % Create ScanList
            app.ScanList = uilistbox(app.ScanListPanel);
            app.ScanList.Items = {'Scans will appear here'};
            app.ScanList.ValueChangedFcn = createCallbackFcn(app, @ScanListValueChanged, true);
            app.ScanList.Position = [19 11 222 260];
            app.ScanList.Value = 'Scans will appear here';

            % Create RelaxometryPanel
            app.RelaxometryPanel = uipanel(app.UIFigure);
            app.RelaxometryPanel.Title = 'Relaxometry';
            app.RelaxometryPanel.Scrollable = 'on';
            app.RelaxometryPanel.Position = [651 136 604 201];

            % Create GenerateRelaxationMapsButton
            app.GenerateRelaxationMapsButton = uibutton(app.RelaxometryPanel, 'push');
            app.GenerateRelaxationMapsButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateRelaxationMapsButtonPushed, true);
            app.GenerateRelaxationMapsButton.Enable = 'off';
            app.GenerateRelaxationMapsButton.Tooltip = {'Produces voxel-wise relaxation parameter maps. See MapRelaxation for details of how the maps are produced'};
            app.GenerateRelaxationMapsButton.Position = [12 146 186 22];
            app.GenerateRelaxationMapsButton.Text = 'Generate Relaxation Maps';

            % Create GenerateROIDispersionCurveButton
            app.GenerateROIDispersionCurveButton = uibutton(app.RelaxometryPanel, 'push');
            app.GenerateROIDispersionCurveButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateROIDispersionCurveButtonPushed, true);
            app.GenerateROIDispersionCurveButton.Enable = 'off';
            app.GenerateROIDispersionCurveButton.Tooltip = {'Produces a ROI based dispersion curve. See T1dispersion for details on how the results are derived.'};
            app.GenerateROIDispersionCurveButton.Position = [12 117 186 22];
            app.GenerateROIDispersionCurveButton.Text = 'Generate ROI Dispersion Curve';

            % Create CheckfitsButton
            app.CheckfitsButton = uibutton(app.RelaxometryPanel, 'state');
            app.CheckfitsButton.Tooltip = {'Displays the relaxation fits when either relaxation maps or ROI dispersion curves are generated'};
            app.CheckfitsButton.Enable = 'off';
            app.CheckfitsButton.Text = 'Check fits';
            app.CheckfitsButton.Position = [215 146 158 22];

            % Create LogarithmicMapsButton
            app.LogarithmicMapsButton = uibutton(app.RelaxometryPanel, 'state');
            app.LogarithmicMapsButton.ValueChangedFcn = createCallbackFcn(app, @LogarithmicMapsButtonValueChanged, true);
            app.LogarithmicMapsButton.Tooltip = {'Visualises relaxation maps either linearly or in log scale. Log scale can be useful for comparing relaxation between different magnetic fields where there is a large change in T1 values'};
            app.LogarithmicMapsButton.Enable = 'off';
            app.LogarithmicMapsButton.Text = 'Logarithmic Maps';
            app.LogarithmicMapsButton.Position = [215 117 158 22];

            % Create DisplayingasT1Button
            app.DisplayingasT1Button = uibutton(app.RelaxometryPanel, 'state');
            app.DisplayingasT1Button.ValueChangedFcn = createCallbackFcn(app, @DisplayingasT1ButtonValueChanged, true);
            app.DisplayingasT1Button.Tooltip = {'Displays dispersion curves or relaxation maps in terms of R1 or T1'};
            app.DisplayingasT1Button.Enable = 'off';
            app.DisplayingasT1Button.Text = 'Displaying as T1';
            app.DisplayingasT1Button.Position = [215 89 158 22];

            % Create ColourAxisUpperLimitEditFieldLabel
            app.ColourAxisUpperLimitEditFieldLabel = uilabel(app.RelaxometryPanel);
            app.ColourAxisUpperLimitEditFieldLabel.HorizontalAlignment = 'right';
            app.ColourAxisUpperLimitEditFieldLabel.Enable = 'off';
            app.ColourAxisUpperLimitEditFieldLabel.Position = [394 149 132 22];
            app.ColourAxisUpperLimitEditFieldLabel.Text = 'Colour Axis Upper Limit';

            % Create ColourAxisUpperLimitEditField
            app.ColourAxisUpperLimitEditField = uieditfield(app.RelaxometryPanel, 'numeric');
            app.ColourAxisUpperLimitEditField.Limits = [0 Inf];
            app.ColourAxisUpperLimitEditField.ValueChangedFcn = createCallbackFcn(app, @ColourAxisUpperLimitEditFieldValueChanged, true);
            app.ColourAxisUpperLimitEditField.Enable = 'off';
            app.ColourAxisUpperLimitEditField.Tooltip = {'Upper limit of displayed T1/R1 values in relaxation maps'};
            app.ColourAxisUpperLimitEditField.Position = [540 149 59 22];
            app.ColourAxisUpperLimitEditField.Value = 1000;

            % Create ColourAxisLowerLimitEditFieldLabel
            app.ColourAxisLowerLimitEditFieldLabel = uilabel(app.RelaxometryPanel);
            app.ColourAxisLowerLimitEditFieldLabel.HorizontalAlignment = 'right';
            app.ColourAxisLowerLimitEditFieldLabel.Enable = 'off';
            app.ColourAxisLowerLimitEditFieldLabel.Position = [394 120 132 22];
            app.ColourAxisLowerLimitEditFieldLabel.Text = 'Colour Axis Lower Limit';

            % Create ColourAxisLowerLimitEditField
            app.ColourAxisLowerLimitEditField = uieditfield(app.RelaxometryPanel, 'numeric');
            app.ColourAxisLowerLimitEditField.Limits = [0 Inf];
            app.ColourAxisLowerLimitEditField.ValueChangedFcn = createCallbackFcn(app, @ColourAxisLowerLimitEditFieldValueChanged, true);
            app.ColourAxisLowerLimitEditField.Enable = 'off';
            app.ColourAxisLowerLimitEditField.Tooltip = {'Lower limit of displayed T1/R1 values in relaxation maps'};
            app.ColourAxisLowerLimitEditField.Position = [540 120 59 22];

            % Create Panel_2
            app.Panel_2 = uipanel(app.UIFigure);
            app.Panel_2.Scrollable = 'on';
            app.Panel_2.Position = [995 655 260 72];

            % Create LoadDataSetButton
            app.LoadDataSetButton = uibutton(app.Panel_2, 'push');
            app.LoadDataSetButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataSetButtonPushed, true);
            app.LoadDataSetButton.Tooltip = {'Open a FFC data set. '};
            app.LoadDataSetButton.Position = [7 38 139 22];
            app.LoadDataSetButton.Text = 'Load Data Set';

            % Create ReprocessingButton
            app.ReprocessingButton = uibutton(app.Panel_2, 'state');
            app.ReprocessingButton.Tooltip = {'Forces the complete reconstruction of the data set. Otherwise the software will look for preprocessed data and open that instead. Can be time consuming for large data sets.'};
            app.ReprocessingButton.Text = 'Reprocessing';
            app.ReprocessingButton.Position = [149 9 101 22];

            % Create BackgroundselectionButton
            app.BackgroundselectionButton = uibutton(app.Panel_2, 'state');
            app.BackgroundselectionButton.Tooltip = {'Allows the user to select the image background to aid in image reconstruction. Usually not needed but may be useful if significant ghosting remains.'};
            app.BackgroundselectionButton.Text = 'Background selection';
            app.BackgroundselectionButton.Position = [7 9 139 22];

            % Create ExportOptionsPanel
            app.ExportOptionsPanel = uipanel(app.UIFigure);
            app.ExportOptionsPanel.AutoResizeChildren = 'off';
            app.ExportOptionsPanel.Title = 'Export Options';
            app.ExportOptionsPanel.Position = [822 10 433 109];

            % Create ExportscantoworkspaceButton
            app.ExportscantoworkspaceButton = uibutton(app.ExportOptionsPanel, 'push');
            app.ExportscantoworkspaceButton.ButtonPushedFcn = createCallbackFcn(app, @ExportscantoworkspaceButtonPushed, true);
            app.ExportscantoworkspaceButton.Enable = 'off';
            app.ExportscantoworkspaceButton.Tooltip = {'Export the scan object to the MATLAB workspace'};
            app.ExportscantoworkspaceButton.Position = [9 63 145 22];
            app.ExportscantoworkspaceButton.Text = 'Export scan to workspace';

            % Create ExportallimagesButton
            app.ExportallimagesButton = uibutton(app.ExportOptionsPanel, 'push');
            app.ExportallimagesButton.ButtonPushedFcn = createCallbackFcn(app, @ExportallimagesButtonPushed, true);
            app.ExportallimagesButton.Enable = 'off';
            app.ExportallimagesButton.Tooltip = {'Export the current data set (magnitude, phase or relaxation map) as a series of .tif images.'};
            app.ExportallimagesButton.Position = [9 34 145 22];
            app.ExportallimagesButton.Text = 'Export all images';

            % Create ExportstackedimagesButton
            app.ExportstackedimagesButton = uibutton(app.ExportOptionsPanel, 'push');
            app.ExportstackedimagesButton.ButtonPushedFcn = createCallbackFcn(app, @ExportstackedimagesButtonPushed, true);
            app.ExportstackedimagesButton.Enable = 'off';
            app.ExportstackedimagesButton.Tooltip = {'Export the current data set (magnitude, phase or relaxation map) as a stacked .tif image. For magnitude and phase images a stacked image of evolution times will be generated for each evolution field. For relaxation maps a single stack containing all relaxation maps will be generated.'};
            app.ExportstackedimagesButton.Position = [9 7 145 22];
            app.ExportstackedimagesButton.Text = 'Export stacked images';

            % Create ExportscanparameterstoworkspaceButton
            app.ExportscanparameterstoworkspaceButton = uibutton(app.ExportOptionsPanel, 'push');
            app.ExportscanparameterstoworkspaceButton.ButtonPushedFcn = createCallbackFcn(app, @ExportscanparameterstoworkspaceButtonPushed, true);
            app.ExportscanparameterstoworkspaceButton.Enable = 'off';
            app.ExportscanparameterstoworkspaceButton.Tooltip = {'Write the current scan parameters in the workspace'};
            app.ExportscanparameterstoworkspaceButton.Position = [160 63 218 22];
            app.ExportscanparameterstoworkspaceButton.Text = 'Export scan parameters to workspace';

            % Create CreateFigurePanel
            app.CreateFigurePanel = uipanel(app.ExportOptionsPanel);
            app.CreateFigurePanel.AutoResizeChildren = 'off';
            app.CreateFigurePanel.Title = 'Create Figure';
            app.CreateFigurePanel.Scrollable = 'on';
            app.CreateFigurePanel.Position = [160 4 218 51];

            % Create kspaceButton
            app.kspaceButton = uibutton(app.CreateFigurePanel, 'push');
            app.kspaceButton.ButtonPushedFcn = createCallbackFcn(app, @kspaceButtonPushed, true);
            app.kspaceButton.Enable = 'off';
            app.kspaceButton.Tooltip = {'Produce a new figure containing the currently displayed k-space matrix.'};
            app.kspaceButton.Position = [5 2 100 22];
            app.kspaceButton.Text = 'k-space';

            % Create ImageButton
            app.ImageButton = uibutton(app.CreateFigurePanel, 'push');
            app.ImageButton.ButtonPushedFcn = createCallbackFcn(app, @ImageButtonPushed, true);
            app.ImageButton.Enable = 'off';
            app.ImageButton.Tooltip = {'Produce a new figure containing the currently displayed image.'};
            app.ImageButton.Position = [109 2 100 22];
            app.ImageButton.Text = 'Image';

            % Create MultichannelDropDown
            app.MultichannelDropDown = uilistbox(app.UIFigure);
            app.MultichannelDropDown.Items = {};
            app.MultichannelDropDown.Multiselect = 'on';
            app.MultichannelDropDown.ValueChangedFcn = createCallbackFcn(app, @MultichannelDropDownValueChanged2, true);
            app.MultichannelDropDown.Enable = 'off';
            app.MultichannelDropDown.Tooltip = {'Multichannel reconstruction options. Currently supports sum of squares and channel averaging. Individual channels can also be inspected.'};
            app.MultichannelDropDown.BackgroundColor = [0.96 0.96 0.96];
            app.MultichannelDropDown.Position = [28 358 63 369];
            app.MultichannelDropDown.Value = {};

            % Create QualityAssurancePanel
            app.QualityAssurancePanel = uipanel(app.UIFigure);
            app.QualityAssurancePanel.Title = 'Quality Assurance';
            app.QualityAssurancePanel.Position = [218 129 417 71];

            % Create GeometricAccButton
            app.GeometricAccButton = uibutton(app.QualityAssurancePanel, 'push');
            app.GeometricAccButton.ButtonPushedFcn = createCallbackFcn(app, @GeometricAccButtonPushed, true);
            app.GeometricAccButton.Position = [9 11 93 23];
            app.GeometricAccButton.Text = 'Geometric Acc';

            % Create UniformityButton
            app.UniformityButton = uibutton(app.QualityAssurancePanel, 'push');
            app.UniformityButton.ButtonPushedFcn = createCallbackFcn(app, @UniformityButtonPushed, true);
            app.UniformityButton.Position = [114 11 85 23];
            app.UniformityButton.Text = 'Uniformity';

            % Create SNRButton
            app.SNRButton = uibutton(app.QualityAssurancePanel, 'push');
            app.SNRButton.ButtonPushedFcn = createCallbackFcn(app, @SNRButtonPushed, true);
            app.SNRButton.Enable = 'off';
            app.SNRButton.Tooltip = {'Attempts to determine the overall image SNR. This method is better suited to QA style measurements with phantoms rather than in-vivo experiments.'};
            app.SNRButton.Position = [211 11 58 23];
            app.SNRButton.Text = 'SNR';

            % Create THKButton
            app.THKButton = uibutton(app.QualityAssurancePanel, 'push');
            app.THKButton.ButtonPushedFcn = createCallbackFcn(app, @THKButtonPushed, true);
            app.THKButton.Position = [280 11 59 23];
            app.THKButton.Text = 'THK';

            % Create AllButton
            app.AllButton = uibutton(app.QualityAssurancePanel, 'push');
            app.AllButton.ButtonPushedFcn = createCallbackFcn(app, @AllButtonPushed, true);
            app.AllButton.Position = [350 11 35 23];
            app.AllButton.Text = 'All';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = FFCProcessing

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            closeFilterDialogs(app);
            if ~isempty(app.cviDialog) && isgraphics(app.cviDialog), delete(app.cviDialog); end
            if ~isempty(app.JobRunner), delete(app.JobRunner); end
            delete(app.UIFigure)
        end
    end
end
