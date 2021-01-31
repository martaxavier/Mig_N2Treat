function [wpli, f_vector, n_wins, n_fft, msgs] = ...
    wpli(X, Y, fs, max_freq_res, overlap, is_symmetric)

%   Computes the connectivity spectrum between signals in X and Y
%   using the weighted phase lag index measure 
%
%   INPUTS:
%
%           X, Y            - Input signals (#chans x #pnts)
%           fs              - Sampling frequency of the input data, X
%           max_freq_res    - Maximum frequency resolution for TF-decomp. 
%           overlap         - % overlap between windows 
%           is_symmetric    - If resulting connectivity matrices are
%                             symmetric
%
%   OUTPUTS:
%
%           wpli            - Connectivity spectrum of the signal 
%           f_vector        - Vector of frequencies of TF-decomp.
%           n_wins          - Number of windows of windowing procedure
%           n_fft           - Number of bins of FFT used for TF-decomp. 
% 

%------------------------------------------------------------
% Initializations
%------------------------------------------------------------

msgs = [];
Cxy = [];
f_vector = [];

% Signal properties
n_X = size(X, 1); 
n_Y = size(Y, 1);
n_pnts = size(X, 2);

% Minimum number of windows
min_win_error = 2;
min_win_warning = 5;

%------------------------------------------------------------
% Frequency Resolution 
%------------------------------------------------------------

% Convert maximum frequency resolution to time length
n_fft = 2^nextpow2( round(fs / max_freq_res) );

% Number of segments
n_overlap = floor(overlap * n_fft);
n_wins = floor((n_pnts - n_overlap) / (n_fft - n_overlap));

% ERROR: Not enough time points
if (n_pnts < n_fft) || (n_wins < min_win_error)
    min_times = n_fft + (n_fft - n_overlap) * (min_win_error - 1);
    msgs = sprintf([strcat('Input signals are too short (%d samples) for', ...
        ' the requested frequency resolution (%1.2fHz).\n') ...
        strcat('Minimum length for this resolution: %1.3f', ...
        ' seconds (%d samples).')], n_pnts, max_freq_res, ...
        min_times/fs, min_times);
    return;
    
% WARNING: Maybe not enough time points
elseif (n_wins < min_win_warning)
    min_times = n_fft + (n_fft - n_overlap) * (min_win_warning - 1);
    msgs = sprintf([strcat('Input signals may be too short (%d samples) for', ...
        ' the requested frequency resolution (%1.2fHz).\n') ...
        strcat('Recommended length for this resolution: %1.3f seconds', ...
        ' (%d samples).')], n_pnts, max_freq_res, min_times/fs, min_times);
end

% Output vector of frequencies
f_vector = fs/2 * linspace(0, 1, n_fft/2 + 1)';
f_vector(end) = [];

%------------------------------------------------------------
% Windowing
%------------------------------------------------------------

% Segment indices - discard final timepoints
i_start = (n_fft - n_overlap) * (0 : (n_wins - 1)) + 1;
i_stop  = i_start + (n_fft-1);
if (i_stop(end) > n_pnts)
    i_start(end) = [];
    i_stop(end) = [];
    n_wins = n_wins - 1;
end

i_win = [i_start; i_stop];

% frequency smoother (represented as time-domain multiplication)
smoother = bst_window('parzen', n_fft) .* bst_window('tukey', n_fft, 0.1);
smoother = smoother / sqrt(sum(smoother.^2));

%------------------------------------------------------------
% Version 2 - symmetric
%------------------------------------------------------------

if is_symmetric
    
    % Indices for the multiplication
    [iY,iX] = meshgrid(1:n_X,1:n_Y);
    
    % Find the values above the diagonal
    ind_sym = find(iX <= iY);
    
    % Cross-spectrum
    Cxy_win = zeros(length(ind_sym), length(f_vector));
    Cxy = zeros(length(ind_sym), length(f_vector), n_wins);
    
    for i = 1 : n_wins
        
        % Get time indices for this segment
        i_time = i_win(1,i) : i_win(2,i);
        
        % frequency domain spectrum after smoothing and tapering
        fourierX = fft(bsxfun(@times, X(:,i_time), smoother'), n_fft, 2);
        fourierY = conj(fft(bsxfun(@times, Y(:,i_time), smoother'), n_fft, 2));
        
        % Calculate for each frequency: fourierX * fourierY'
        Cxy_win = fourierX(iX(ind_sym),1:(n_fft/2)) ...
            .* fourierY(iY(ind_sym),1:(n_fft/2));
        Cxy(:, :, i) = Cxy_win;
        
    end
           
    imagsum      = sum(imag(Cxy),3);
    imagsumW     = sum(abs(imag(Cxy)),3);
    debiasfactor = sum(imag(Cxy).^2,3);
    wpli = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);   
    
end

    % Make sure that there are
    % no residual imaginary parts
    % due to numerical errors
    if ~isreal(wpli)
        wpli = abs(wpli);
    end

end
