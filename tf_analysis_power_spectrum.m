function [power,f_vector] = tf_analysis_power_spectrum(data, ...
    f_range, n_freq, tf_method, wav_seconds, win_seconds, ...
    n_wins, fs_data, fs_power)

%   Performs time-frequency analysis of the signals  in 'data'
%   Extracts the power-spectrum of each signal in 'data', using 
%   the method 'tf_method'
%
%   INPUTS:
%
%           data        - Input signals 
%           f_range     - Frequency range for the TF-decomposition 
%           n_freq      - Number of frequency bins for the TF-decomposition 
%           tf_method   - Method to be used in the TF-decomposition 
%           wav_seconds - Size of the Morlet Wavelet kernel (seconds)
%           win_seconds - Size of the Welch sliding window (seconds)
%           fs_data     - Sampling frequency of the input data, 'data'
%           fs_power    - Sampling frequency of the output data, 'power'
% 
%   OUTPUTS:
%
%           power       - Output power time-series of each input signal 
%           f_vector    - Frequency vector used in thr tf-decomposition 
% 

% Ensure that data is a real 2D matrix
if ~ismatrix(data) || ~isreal(data)
    error('data is not a real 2D matrix');
end

% Ensure that the first dimension 
% spans time-points 
if size(data,1) < size(data,2)
    data = squeeze(data');
end

% Number of time-points
n_pnts = size(data, 1);

% Obtain number of time points of the output 
n_pnts_power = round((n_pnts - 1)*(fs_power/fs_data) + 1);

switch tf_method
    
    % === Morlet Wavelet ===
    case 'wavelet'
        
        % Extract the tf power spectrum of the current
        % channel, using Morlet wavelet decomposition 
        [~ , f_vector, ~, power] = tf_power_spectrum_wavelet(data, ...
            f_range, n_freq, wav_seconds, fs_data, 1);
        
    % === Welch's Method ===       
    case 'welch'

        % Extract the tf power spectrum of the current 
        % channel, using the Welch method 
        % (the cross-spectrum of one signal with itself
        % is its power spectral density)
        [f_vector, power] =  ...
            tf_cross_spectrum_welch(data, data, f_range, n_freq, ...
            n_wins, win_seconds, fs_data, fs_power, 0);

        % Padd the tf power-spectrum to match the begining and 
        % end of the input data 
        padd_size = (n_pnts_power - size(power, 1)) / 2;
        power = padarray(power, [padd_size 0 0], ...
            'replicate','both');  

end

end