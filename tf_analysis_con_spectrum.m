function [Cxy, p_values, f_vector] = tf_analysis_con_spectrum(X, ...
    f_range, n_freq, win_seconds, fs_X, fs_Cxy, metric, stat_filt_method)

%   Performs time-frequency analysis of the signals  in X
%   Computes the connectivity spectrum between each pair
%   of signals in X 
%
%   INPUTS:
%
%           X           - Input signals (#chans x #pnts)
%           f_range     - Frequency range for the TF-decomposition 
%           n_freq      - Number of frequency bins for the TF-decomposition 
%           win_seconds - Size of the Welch sliding window (seconds)
%           fs_X        - Sampling frequency of the input data, X
%           fs_Cxy      - Sampling frequency of the output data, Cxy
%           metric      - Connectivity metric 
% 
%   OUTPUTS:
%
%           Cxy         - Connectivity spectrum of the signal 
% 

%------------------------------------------------------------
% Read and check input data 
%------------------------------------------------------------

% Ensure that the first 
% dimension spans channels 
if size(X,1) > size(X,2)
    X = squeeze(X');
end
    
% Read size of input data 
n_pnts_X = size(X,2);

%------------------------------------------------------------
% Parameters of the TF decomposition 
%------------------------------------------------------------

% Parameters for TF decomposition 
win_step = fs_X/fs_Cxy;                         % window step (samples)
win_samples = win_seconds*fs_X + 1;             % window size (samples)

% Define frequency resolution 
f_max = f_range(2);                             % maximum frequency 
f_min = f_range(1);                             % minimum frequency

f_vector = linspace(f_min, f_max, n_freq);      % vector of frequencies
max_freq_res = max(diff(f_vector));             % maximum frequency 
                                                % resolution of the vector
                                                % of frequencies 
                                                
%------------------------------------------------------------
% Use sliding window to compute dynamic connectivity 
%------------------------------------------------------------


% Obtain number of time-points of the output 
n_pnts_Cxy = round((n_pnts_X - 1)*(fs_Cxy/fs_X) + 1);

% Discard the first and last time-points 
start_X = 1 + ((win_samples - 1) / 2);
stop_X = n_pnts_X - ((win_samples - 1) / 2);

% Vector containing the center of each window,  
% in the fs of the input data  
win_center_X = round(start_X : win_step : stop_X);

% Compute connectivity for the first window center
% to pre-allocate the connectivity and p-value matrices 
win = win_center_X(1) - (win_samples - 1) / 2 ...
    : win_center_X(1) + (win_samples - 1) / 2;
X_win = X(:, win);
[Cxy_win, p_values_win, f_vector] = ...
    get_connectivity(X_win, fs_X, max_freq_res, metric);

Cxy = zeros([length(win_center_X), size(Cxy_win)]);
Cxy(1, :, :) = Cxy_win; 

if strcmp(metric,'icoh') && ...
        strcmp(stat_filt_method, 'analytical')
    p_values = Cxy;
    p_values(1, :, :) = p_values_win;
else
    p_values = [];
end

center_Cxy = 2; 

% Go through each window 
for center_X = win_center_X(2 : end)

    % Define the current window, centered on 'center_Xs'
    win = center_X - (win_samples - 1) / 2 ...
        : center_X + (win_samples - 1) / 2;
    X_win = X(:, win);

    % Compute the specified connectivity metric for the
    % current window and update vector of frequencies
    [Cxy_win, p_values_win, f_vector] = get_connectivity(X_win, ...
        fs_X, max_freq_res, metric);
    Cxy(center_Cxy, :, :) = Cxy_win;
    
    if strcmp(metric, 'icoh')&& ...
            strcmp(stat_filt_method, 'analytical')
        p_values(center_Cxy, :, :) = p_values_win; 
    end
    
    % Update the time-point for which the 
    % power spectral density is obtained 
    center_Cxy = center_Cxy + 1;

end

% Padd the connectivity to match the begining and 
% spectrum to match the end of the input data 
padd_size = (n_pnts_Cxy - length(win_center_X)) / 2;
Cxy = padarray(Cxy, [padd_size 0], ...
    'replicate','both');

if strcmp(metric,'icoh') && ...
        strcmp(stat_filt_method, 'analytical')
    p_values = padarray(p_values, ...
        [padd_size 0], 'replicate','both');
end

end

    function [Cxy, p_values, f_vector] = ...
        get_connectivity(X, fs_X, max_freq_res, metric)
    
    p_values = [];
    
    switch metric
        
        % Imaginary part of coherence
        case 'icoh'
            
            [Cxy, p_values, f_vector, ~, ~, ~] = ...
                bst_cohn(X, X, fs_X, max_freq_res, 0.5, ...
                'icohere', 1, [], 100);
            Cxy = squeeze(Cxy);
        
        % Weighted phase lag index 
        % (debiased) 
        case 'wpli'
            
            [Cxy, f_vector, ~, ~, ~] = ...
                wpli(X, X, fs_X, max_freq_res, 0.5, 1);
            Cxy = squeeze(Cxy);
            
    end
        
    end