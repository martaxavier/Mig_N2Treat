function [cspec, f_vector] = tf_analysis_cross_spectrum(data, ...
    f_range, n_freq, n_wins_welch, win_seconds, fs_data, fs_cspec)

%   Performs time-frequency analysis of the signals  in 'data'
%   Extracts the power-spectrum of each signal in 'data', using 
%   the method 'tf_method'
%
%   INPUTS:
%
%           data        - Input signals 
%           f_range     - Frequency range for the TF-decomposition 
%           n_freq      - Number of frequency bins for the TF-decomposition 
%           win_seconds - Size of the Welch sliding window (seconds)
%           fs_data     - Sampling frequency of the input data, 'data'
%           fs_cspec    - Sampling frequency of the output data, 'power'
% 
%   OUTPUTS:
%
%           cspec       - Cross-spectrum time-series of each input signal 
%           f_vector    - Frequency vector used in thr tf-decomposition 
% 

% Ensure that data is a real 2D matrix
if ~ismatrix(data) ||  length(size(data)) ~= 2 || ~isreal(data)
    error('Input data is not a real 2D matrix');
end

% Ensure that the first dimension 
% spans channels 
if size(data,1) > size(data,2)
    data = squeeze(data');
end
    
% Read size of input data 
n_chans = size(data,1);
n_pnts = size(data,2);

% Obtain number of time-points of the output 
n_pnts_cspec = round((n_pnts - 1)*(fs_cspec/fs_data) + 1);

% Pre-allocate the output matrix 
cspec = zeros(n_freq, n_pnts_cspec, n_chans, n_chans);

% Go through channels 
for channel1 = 1 : n_chans
    
    % Lower triangular matrix 
    for channel2 = 1 : channel1

        % Extract the tf cross-spectrum of the current 
        % pair of channels, using the Welch method 
        [f_vector, channel_cspec] =  ...
            tf_cross_spectrum_welch(data(channel1,:), ...
            data(channel2,:), f_range, n_freq, n_wins_welch, ...
            win_seconds, fs_data, fs_cspec, 0);

        % Padd the tf cross-spectrum to match the begining and 
        % end of the input data 
        padd_size = (n_pnts_cspec - size(channel_cspec,2)) / 2;
        channel_cspec = padarray(channel_cspec, [0 padd_size], ...
            'replicate','both');
        
        % NOTE dof will be exactly the same for every 
        % element of the matrix; we should define dof outside 
        % and we should define the parameters outside as well 
        cspec(:,:,channel1,channel2) = channel_cspec;
            
    end
    
end

end