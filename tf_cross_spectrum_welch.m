function [f_vector, Cxy] = tf_cross_spectrum_welch(X, Y, ...
    f_range, n_freq, n_wins, win_seconds, fs_X, fs_Cxy, plot_TF)
%
%   function [f_vector, power] = tf_cross_spectrum_welch(X,
%   f_range n_freq, fs_data, fs_power, plot_TF) computes the 
%   cross-spectrum of the signals X and Y, using the Welch method 
%
%   INPUTS:
%
%           X, Y        - Input signals 
%           f_range     - Frequency range for the TF-decomposition 
%           n_freq      - Number of frequency bins for the TF-decomposition
%           n_wins      - Number of windows for Welch TF-decomposition 
%           win_seconds - Sliding window size(in seconds)
%           fs_X        - Sampling frequency of the input data, 'X'/'Y'
%           fs_Cxy      - Sampling frequency of the output data, 'Cxy' 
%           plot_TF     - 1 to plot, 0 to don't 
%
%   OUTPUTS:
%
%           f_vector    - Frequency vector used in thr tf-decomposition
%           Cxy         - Cross-spectrum time-series of each input signal 

%------------------------------------------------------------
% Prepare input for TF decomposition 
%------------------------------------------------------------

% Parameters for TF decomposition 
win_step = fs_X/fs_Cxy;                         % window step (samples)
win_samples = win_seconds*fs_X + 1;             % window size (samples)

% Build vector of frequencies 
f_max = f_range(2);                             % maximum frequency 
f_min = f_range(1);                             % minimum freuency

f_vector = linspace(f_min,f_max,n_freq);        % vector of frequencies

max_freq_res = max(diff(f_vector));             % maximum frequency 
                                                % resolution of the vector
                                                % of frequencies 
                                                
n_fft = 2*floor(win_samples / n_wins);          % number of FFT samples
                                                % the factor 2 is because
                                                % there is 50% overlap

if n_fft < round(fs_X / max_freq_res) 
    warning(strcat('The number of FFT samples as', ...
        ' specified by the number of Welch windows', ...
        ' is lower than number of FFT samples', ...
        ' necessary for the specified frequency', ...
        ' resolution. MATLAB will perform zero', ...
        ' padding before computing the FFT'));
end           

if n_wins < 5
    warning('Number of Welch windows is too short');
end

%-------------------------------------------------------------------------%
% Perform TF decomposition using the Welch method 
%-------------------------------------------------------------------------%

% Number of time-points 
% of the input data 
n_pnts_X = length(X);

% Discard the first and last time-points 
start_X = 1 + ((win_samples - 1) / 2);
stop_X = n_pnts_X - ((win_samples - 1) / 2);

% Vector containing the center of each window,  
% in the fs of the input data  
win_center_X = round(start_X : win_step : stop_X);

% Number of time-points of the output data 
n_pnts_Cxy = length(win_center_X);

% Center of the first window, in
% the fs of the output data 
center_Cxy = 1; 

% Pre-allocate matrix containing the power
% spectral densitity througout time 
Cxy = zeros(n_freq,n_pnts_Cxy);

% Go through each window 
for center_X = win_center_X

    % Define the current window, centered on 'center_Xs'
    win = center_X - (win_samples - 1) / 2 ...
        : center_X + (win_samples - 1) / 2;
    X_win = X(win);
    Y_win = Y(win);
    
    % Obtain the power spectral density for the current window 
    % (the cross-power spectral density - cpsd - of one signal  
    % with itself is its power spectral density)
    % Tapering - Hamming window of 250 ms; Window overlap - 50% 
    Cxy(:,center_Cxy) = cpsd(X_win, Y_win, ...
        hamming(n_fft), [], f_vector, fs_X);
    
    % Update the time-point for which the 
    % power spectral density is obtained 
    center_Cxy = center_Cxy + 1;

end

% Plot power spectrum
if plot_TF
    
    % Colorscale
    Dtf_f = log(abs(Cxy) + .001); 
    
    figure('Name', 'TF')

    imagesc((1:size(Cxy,2)) ./ fs_Cxy, ...
        1:length(f_vector), squeeze(Dtf_f));
    
%     Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) ...
%         - log(f_vector(1))) ./ log(1 + alpha) + 2;
%     hold on; plot([1 n_pnts_Cxy],[Fplot', Fplot'],'k');
%     hold off
    
    set(gca,'YTick',Fplot);
    set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90],'FontSize',12)
    ylabel('Frequency (Hz)','FontSize', 26);
    xlabel('Time (s)','FontSize',26);
    colorbar;

end