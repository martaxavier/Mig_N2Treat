function [f_vector, Cxy] = tf_cross_spectrum_welch(X,Y, ...
    f_range, n_freq, win_seconds, fs_X, fs_power, plot_TF)
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
%           win_seconds - Sliding window size(in seconds)
%           fs_X        - Sampling frequency of the input data, 'X'
%           fs_cspec    - Sampling frequency of the output data, 'cspec' 
%           plot_TF     - 1 to plot, 0 to don't 
%
%   OUTPUTS:
%
%           cspec       - Cross-spectrum time-series of each input signal 
%           f_vector    - Frequency vector used in thr tf-decomposition
%           dof         - 2.8 * win_samples * (1.28/(nFFT+1))
%                         [Welch 1976]

%------------------------------------------------------------
% Prepare input for TF decomposition 
%------------------------------------------------------------

% Parameters for TF decomposition 
win_step = round(fs_X/fs_power);                % window step (samples)
win_samples = win_seconds*fs_X + 1;             % window size (samples)

% Build vector of frequencies 
f_max = f_range(2);                             % maximum frequency 
f_min = f_range(1);                             % minimum freuency

alpha = (f_max/f_min)^(1/n_freq) - 1; 
f_vector = ((1+alpha).^(1:n_freq)) * f_min;     % vector of frequencies 
                                                % (more frequency bins for
                                                % the lower frequencies)

%-------------------------------------------------------------------------%
% Perform TF decomposition using the Welch method 
%-------------------------------------------------------------------------%

% Vector of time-points, in the output fs 
n_pnts_X = length(X);
%t_vector_power = 0 : 1/fs_power : (n_pnts_X -1)/fs_X;
%n_pnts_power = length(t_vector_power);

% Discard the first and last time-points 
%win_offset_X = (win_samples - 1) / 2;
start_X = 1 + ((win_samples - 1) / 2);
stop_X = n_pnts_X - ((win_samples - 1) / 2);

% Vector containing the center of each window,  
% in the input fs 
win_center_X = round(start_X : win_step : stop_X);
n_pnts_power = length(win_center_X);

% Center of the first window, in the output fs
center_power = 1; 
%center_power = 1 + (n_pnts_power - length(win_center_X)) / 2;

% Pre-allocate matrix containing the power
% spectral densitity througout time 
Cxy = zeros(n_freq,n_pnts_power);

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
    Cxy(:,center_power) = cpsd(X_win, Y_win, ...
        hamming(250), [], f_vector, fs_X);
    
    % Update the time-point for which the 
    % power spectral density is obtained 
    center_power = center_power + 1;

end

% Obtain the window offset, in the output fs 
%win_offset_power = (win_offset_X - 1) * (fs_power/fs_X) + 1;

% Padd the output signal, outside the samples 
% for which the power spectral density was obtained 
%power = Cxy(:, round(1 + win_offset_power : ...
%    end - win_offset_power));
%power = padarray(power, [0 win_offset_power], ...
%    'replicate','both');

% Plot power spectrum
if plot_TF
    
    % Colorscale
    Dtf_f = log(abs(power) + .001); 
    
    figure('Name', 'TF')

    imagesc((1:size(power,2)) ./ fs_X, ...
        1:length(f_vector), squeeze(Dtf_f));
    
    Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) ...
        - log(f_vector(1))) ./ log(1 + alpha) + 2;
    hold on; plot([1 n_pnts_power],[Fplot', Fplot'],'k');
    hold off
    
    set(gca,'YTick',Fplot);
    set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90],'FontSize',12)
    ylabel('Frequency (Hz)','FontSize', 26);
    xlabel('Time (s)','FontSize',26);
    colorbar;

end