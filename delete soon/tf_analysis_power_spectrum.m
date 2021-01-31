function [power,f_vector] = tf_analysis_power_spectrum(data, ...
    f_range, n_freq, tf_method, wav_seconds, win_seconds, ...
    fs_data, fs_power)

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
% spans channels 
if size(data,1) > size(data,2)
    data = squeeze(data');
end
    
% Read size of input data 
n_chans = size(data,1);
n_pnts = size(data,2);

% Obtain number of time points of the output 
n_pnts_power = round((n_pnts - 1)*(fs_power/fs_data) + 1);

% Pre-allocate the output matrix 
switch tf_method
    case 'wavelet'
        power = zeros(n_freq,n_pnts,n_chans);
    case 'welch'
        power = zeros(n_freq,n_pnts_power,n_chans);
end

% Go through channels 
for channel = 1 : n_chans
    
    switch tf_method
        
        % === Morlet Wavelet === 
        case 'wavelet'
            
            % Extract the tf power spectrum of the current
            % channel, using Morlet wavelet decomposition 
            [~ , f_vector, ~, channel_pow] = ...
                tf_power_spectrum_wavelet(data(channel,:), ...
                f_range, n_freq, wav_seconds, ...
                fs_data, 1);
            
        % === Welch ===   
        case 'welch'
            
            % Extract the tf power spectrum of the current 
            % channel, using the Welch method 
            % (the cross-spectrum of one signal with itself
            % is its power spectral density)
            [f_vector, channel_pow] =  ...
                tf_cross_spectrum_welch(data(channel,:), ...
                data(channel,:),f_range, n_freq, ...
                win_seconds, fs_data, fs_power, 0);
            
            % Padd the tf power-spectrum to match the begining and 
            % end of the input data 
            padd_size = (n_pnts_power - size(channel_pow,2)) / 2;
            channel_pow = padarray(channel_pow, [0 padd_size], ...
                'replicate','both');            
            
    end
    
    power(:,:,channel) = channel_pow;
    
end

end