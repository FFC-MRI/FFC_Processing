function correctedkspace = correct_phase(kspace, backgroundtest, n_receivers)
%CORRECT_PHASE  Estimate ky-dependent phase error from one coil and apply to all coils.
%
% correctedkspace = correct_phase(kspace, backgroundtest, n_receivers)
%
% kspace expected size: [Nx Ny ... Nc] (channels in last dim)
% n_receivers: number of coils (Nc)
%
% Linear ky ordering assumed.

    dims = size(kspace);
    if numel(dims) < 3
        correctedkspace = kspace;
        return;
    end

    Nx = dims(1);
    Ny = dims(2);
    Nc = dims(end);

    if nargin < 3 || isempty(n_receivers)
        n_receivers = Nc;
    end
    n_receivers = min(n_receivers, Nc);

    % --- Pick coil used for estimating the correction ---
    refCoil = n_receivers; % matches your original intent
    kRef = kspace;
    % Extract only that coil, keep all other dims
    idx = repmat({':'}, 1, ndims(kspace));
    idx{ndims(kspace)} = refCoil;
    kRef = kspace(idx{:});                 % [Nx Ny ...] (no coil dim)

    % Collapse pages (everything beyond Nx,Ny) into 3rd dim: [Nx Ny Npages]
    kRefPages = reshape(kRef, Nx, Ny, []);

    % --- Optional user-defined background mask on one representative image ---
    if backgroundtest == 1
        hh = figure;
        imagesc(abs(fft2c(kRefPages(:,:,1)))); axis off square; colormap('gray');
        background = roipoly;
        close(hh);
    else
        background = [];
    end

    % Estimate phase per page (iterative_images_correction_v7 supports stacks)
    % Output ph will be [1 Ny Npages]
    [~, ph] = iterative_images_correction_v9(kRefPages, 1, 50, 1e-3, background);

    % Ensure ph is [1 Ny Npages]
    if ismatrix(ph)
        ph = reshape(ph, 1, Ny, []);
    end

    % Build phase correction array shaped for implicit expansion:
    % We want: [Nx Ny Npages] so exp(-1i*...) can multiply pages.
    % ph: [1 Ny Npages] -> expand along Nx
    phasePages = reshape(ph, 1, Ny, size(kRefPages,3));
    phasePages = repmat(phasePages, Nx, 1, 1);   % [Nx Ny Npages]

    % Apply to ALL coils at once:
    % Reshape full kspace to [Nx Ny Npages Nc] then multiply
    kAll = reshape(kspace, Nx, Ny, [], Nc);      % [Nx Ny Npages Nc]
    kAll = kAll .* exp(-1i * phasePages);        % implicit expansion over Nc

    correctedkspace = reshape(kAll, dims);
     % correctedkspace = kspace;
end


