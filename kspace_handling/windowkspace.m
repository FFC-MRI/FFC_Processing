function windowed_kspace = windowkspace(kspace, window_size, windowtype)
%WINDOWKSPACE Apply a centred 2D window to k-space.
%
%   windowed_kspace = windowkspace(kspace, window_size, windowtype)
%
%   kspace      : k-space array. First two dimensions are read/phase.
%                 Any trailing dimensions are preserved.
%   window_size : Fractional window size, e.g. 1.0 = full matrix,
%                 0.8 = 80% of matrix.
%   windowtype  : 'None', 'Hann', 'Blackman', 'Hamming',
%                 'Kaiser-Bessel', or 'Tukey'.
%
%   The function is robust to odd/even matrix sizes such as 175 x 174.

dims = size(kspace);
n1 = dims(1);
n2 = dims(2);

kspace_reshaped = reshape(kspace, n1, n2, []);

% Ensure sensible window size
if nargin < 2 || isempty(window_size)
    window_size = 1;
end

if nargin < 3 || isempty(windowtype)
    windowtype = 'None';
end

window_size = max(window_size, 0);

% Calculate requested window dimensions.
% Do NOT force even sizes. Odd-sized matrices need odd-sized windows.
windowdim1 = max(1, round(window_size * n1));
windowdim2 = max(1, round(window_size * n2));

% Choose window function
switch windowtype
    case 'None'
        matlab_window = 'rectwin';
        opts = [];
    case 'Hann'
        matlab_window = 'hann';
        opts = [];
    case 'Blackman'
        matlab_window = 'blackman';
        opts = [];
    case 'Hamming'
        matlab_window = 'hamming';
        opts = [];
    case 'Kaiser-Bessel'
        matlab_window = 'kaiser';
        opts = 3;
    case 'Tukey'
        matlab_window = 'tukeywin';
        opts = 0.5;
    otherwise
        matlab_window = 'rectwin';
        opts = [];
end

% Generate 1D windows
if isempty(opts)
    window1 = window(matlab_window, windowdim1);
    window2 = window(matlab_window, windowdim2);
else
    window1 = window(matlab_window, windowdim1, opts);
    window2 = window(matlab_window, windowdim2, opts);
end

% Outer product to create 2D window
windowfull = window1 * window2.';

% Centre crop or centre pad to exactly match k-space size
windowfull = centre_crop_or_pad_2d(windowfull, n1, n2);

% Apply window, preserving trailing dimensions
windowed_kspace = reshape(kspace_reshaped .* windowfull, dims);

end


function out = centre_crop_or_pad_2d(in, target_n1, target_n2)
%CENTRE_CROP_OR_PAD_2D Centre crop or zero-pad a 2D array to requested size.

[in_n1, in_n2] = size(in);

% First crop if needed
crop_start_1 = floor((in_n1 - min(in_n1, target_n1)) / 2) + 1;
crop_start_2 = floor((in_n2 - min(in_n2, target_n2)) / 2) + 1;

crop_end_1 = crop_start_1 + min(in_n1, target_n1) - 1;
crop_end_2 = crop_start_2 + min(in_n2, target_n2) - 1;

cropped = in(crop_start_1:crop_end_1, crop_start_2:crop_end_2);

% Then pad if needed
out = zeros(target_n1, target_n2, class(in));

[crop_n1, crop_n2] = size(cropped);

insert_start_1 = floor((target_n1 - crop_n1) / 2) + 1;
insert_start_2 = floor((target_n2 - crop_n2) / 2) + 1;

insert_end_1 = insert_start_1 + crop_n1 - 1;
insert_end_2 = insert_start_2 + crop_n2 - 1;

out(insert_start_1:insert_end_1, insert_start_2:insert_end_2) = cropped;

end