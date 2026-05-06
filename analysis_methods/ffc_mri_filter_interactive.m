function filtered_images = ffc_mri_filter_interactive(images)
%FFC_MRI_FILTER_INTERACTIVE Interactive preview/tuning GUI for ffc_mri_filter
%
% Usage:
%   filtered_images = ffc_mri_filter_interactive(images);
%
% images may be [X Y ...]. Extra dimensions are flattened for preview and
% restored on output.

dim = size(images);
temp0 = abs(reshape(images, dim(1), dim(2), []));  % [X Y N]
N = size(temp0, 3);

% Default settings
state.filter_type = 'LLR+Wavelet';
state.frame = 1;

state.patch = 5;
state.step = 3;
state.niter = 2;
state.tau = 0.5;

state.wavelet_level = 3;

state.use_tgv = false;
state.tgv_alpha0 = 0.04;
state.tgv_alpha1 = 0.015;
state.tgv_iter = 25;
state.tgv_maxiter = 150;

state.filtered_preview = [];

% Figure
fig = uifigure('Name','FFC MRI filter tuning', ...
    'Position',[100 100 1200 650]);

gl = uigridlayout(fig, [1 2]);
gl.ColumnWidth = {'1x', 330};

% Image panel
imagePanel = uipanel(gl, 'Title','Preview');
imageLayout = uigridlayout(imagePanel, [2 2]);
imageLayout.RowHeight = {'1x', 40};
imageLayout.ColumnWidth = {'1x','1x'};

ax1 = uiaxes(imageLayout);
title(ax1, 'Original');
axis(ax1, 'image');
ax1.XTick = [];
ax1.YTick = [];

ax2 = uiaxes(imageLayout);
title(ax2, 'Filtered preview');
axis(ax2, 'image');
ax2.XTick = [];
ax2.YTick = [];

statusLabel = uilabel(imageLayout, ...
    'Text','Ready', ...
    'FontWeight','bold');
statusLabel.Layout.Row = 2;
statusLabel.Layout.Column = [1 2];

% Control panel
ctrlPanel = uipanel(gl, 'Title','Filter controls');
ctrl = uigridlayout(ctrlPanel, [18 2]);
ctrl.RowHeight = repmat({30}, 1, 18);
ctrl.ColumnWidth = {120, '1x'};

% Filter dropdown
uilabel(ctrl, 'Text','Filter');
filterDrop = uidropdown(ctrl, ...
    'Items',{'LLR+Wavelet','LLR','Wavelet','LLR+Wavelet+TGV','Total Variation'}, ...
    'Value','LLR+Wavelet');

% Frame slider
uilabel(ctrl, 'Text','Frame');
frameSlider = uislider(ctrl, ...
    'Limits',[1 N], ...
    'Value',1);
frameSlider.MajorTicks = unique(round(linspace(1, N, min(N, 6))));
frameSlider.MinorTicks = [];

% Patch size
uilabel(ctrl, 'Text','Patch');
patchSpin = uispinner(ctrl, ...
    'Limits',[3 32], ...
    'Value',state.patch, ...
    'Step',1);

% Step size
uilabel(ctrl, 'Text','Step');
stepSpin = uispinner(ctrl, ...
    'Limits',[1 16], ...
    'Value',state.step, ...
    'Step',1);

% Iterations
uilabel(ctrl, 'Text','LLR iterations');
niterSpin = uispinner(ctrl, ...
    'Limits',[1 10], ...
    'Value',state.niter, ...
    'Step',1);

% Tau
uilabel(ctrl, 'Text','LLR tau');
tauSlider = uislider(ctrl, ...
    'Limits',[0 5], ...
    'Value',state.tau);
tauSlider.MajorTicks = 0:1:5;

% Wavelet level
uilabel(ctrl, 'Text','Wavelet level');
waveletSpin = uispinner(ctrl, ...
    'Limits',[1 5], ...
    'Value',state.wavelet_level, ...
    'Step',1);

% TGV checkbox
uilabel(ctrl, 'Text','Use TGV');
tgvCheck = uicheckbox(ctrl, ...
    'Value',state.use_tgv);

% TGV alpha0
uilabel(ctrl, 'Text','TGV alpha0');
tgvAlpha0 = uislider(ctrl, ...
    'Limits',[0 0.2], ...
    'Value',state.tgv_alpha0);

% TGV alpha1
uilabel(ctrl, 'Text','TGV alpha1');
tgvAlpha1 = uislider(ctrl, ...
    'Limits',[0 0.1], ...
    'Value',state.tgv_alpha1);

% TGV iterations
uilabel(ctrl, 'Text','TGV iter');
tgvIterSpin = uispinner(ctrl, ...
    'Limits',[1 100], ...
    'Value',state.tgv_iter, ...
    'Step',1);

% Buttons
previewButton = uibutton(ctrl, ...
    'Text','Update preview', ...
    'ButtonPushedFcn',@(~,~) updatePreview());
previewButton.Layout.Column = [1 2];

applyButton = uibutton(ctrl, ...
    'Text','Apply to full stack and close', ...
    'FontWeight','bold', ...
    'ButtonPushedFcn',@(~,~) applyAndClose());
applyButton.Layout.Column = [1 2];

cancelButton = uibutton(ctrl, ...
    'Text','Cancel / return original magnitude', ...
    'ButtonPushedFcn',@(~,~) cancelAndClose());
cancelButton.Layout.Column = [1 2];

% Value labels
tauLabel = uilabel(ctrl, 'Text','');
tauLabel.Layout.Column = [1 2];

tgvLabel = uilabel(ctrl, 'Text','');
tgvLabel.Layout.Column = [1 2];

% callbacks
filterDrop.ValueChangedFcn = @(~,~) updatePreview();
frameSlider.ValueChangedFcn = @(~,~) updatePreview();

patchSpin.ValueChangedFcn = @(~,~) updatePreview();
stepSpin.ValueChangedFcn = @(~,~) updatePreview();
niterSpin.ValueChangedFcn = @(~,~) updatePreview();
tauSlider.ValueChangedFcn = @(~,~) updatePreview();
waveletSpin.ValueChangedFcn = @(~,~) updatePreview();

tgvCheck.ValueChangedFcn = @(~,~) updatePreview();
tgvAlpha0.ValueChangedFcn = @(~,~) updatePreview();
tgvAlpha1.ValueChangedFcn = @(~,~) updatePreview();
tgvIterSpin.ValueChangedFcn = @(~,~) updatePreview();

% Initial display
filtered_images = [];
updatePreview();

% Wait for user to press Apply or Cancel
uiwait(fig);

    function readControls()
        state.filter_type = filterDrop.Value;
        state.frame = round(frameSlider.Value);

        state.patch = patchSpin.Value;
        state.step = stepSpin.Value;
        state.niter = niterSpin.Value;
        state.tau = tauSlider.Value;

        state.wavelet_level = waveletSpin.Value;

        state.use_tgv = tgvCheck.Value;
        state.tgv_alpha0 = tgvAlpha0.Value;
        state.tgv_alpha1 = tgvAlpha1.Value;
        state.tgv_iter = tgvIterSpin.Value;

        % Enforce odd patch size
        if mod(state.patch, 2) == 0
            state.patch = state.patch + 1;
            patchSpin.Value = state.patch;
        end

        tauLabel.Text = sprintf('LLR tau = %.3f', state.tau);
        tgvLabel.Text = sprintf('TGV alpha0 = %.4f, alpha1 = %.4f', ...
            state.tgv_alpha0, state.tgv_alpha1);
    end

    function updatePreview()
        readControls();

        statusLabel.Text = 'Filtering preview...';
        drawnow;

        img = temp0(:,:,state.frame);

        opts = struct();
        opts.patch = state.patch;
        opts.step = state.step;
        opts.niter = state.niter;
        opts.tau = state.tau;
        opts.sigma = [];
        opts.wavelet_level = state.wavelet_level;
        opts.use_tgv = state.use_tgv;
        opts.tgv_alpha0 = state.tgv_alpha0;
        opts.tgv_alpha1 = state.tgv_alpha1;
        opts.tgv_iter = state.tgv_iter;
        opts.tgv_maxiter = state.tgv_maxiter;

        try
            state.filtered_preview = ffc_filter_apply(img, state.filter_type, opts);

            imagesc(ax1, img);
            axis(ax1, 'image');
            colormap(ax1, gray);
            ax1.XTick = [];
            ax1.YTick = [];
            title(ax1, sprintf('Original frame %d/%d', state.frame, N));

            imagesc(ax2, state.filtered_preview);
            axis(ax2, 'image');
            colormap(ax2, gray);
            ax2.XTick = [];
            ax2.YTick = [];
            title(ax2, 'Filtered preview');

            clim = [min(img(:)), max(img(:))];
            if clim(2) > clim(1)
                ax1.CLim = clim;
                ax2.CLim = clim;
            end

            statusLabel.Text = 'Ready';

        catch ME
            statusLabel.Text = ['Error: ' ME.message];
        end

        drawnow;
    end

    function applyAndClose()
        readControls();

        statusLabel.Text = 'Applying filter to full stack...';
        drawnow;

        opts = struct();
        opts.patch = state.patch;
        opts.step = state.step;
        opts.niter = state.niter;
        opts.tau = state.tau;
        opts.sigma = [];
        opts.wavelet_level = state.wavelet_level;
        opts.use_tgv = state.use_tgv;
        opts.tgv_alpha0 = state.tgv_alpha0;
        opts.tgv_alpha1 = state.tgv_alpha1;
        opts.tgv_iter = state.tgv_iter;
        opts.tgv_maxiter = state.tgv_maxiter;

        temp = ffc_filter_apply(temp0, state.filter_type, opts);

        filtered_images = reshape(temp, dim);

        uiresume(fig);
        delete(fig);
    end

    function cancelAndClose()
        filtered_images = reshape(temp0, dim);
        uiresume(fig);
        delete(fig);
    end

end