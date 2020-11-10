function [power_all, freq_vector] = ...
    extract_power(data,f_range,n_freq,tf_method,fs_eeg,fs)

%   Performs time-frequency decomposition on the
%   data to extract power and averages power across
%   the frequency bands defined in the variable bands

%------------------------------------------------------------
% Go through channels 
%------------------------------------------------------------
n_chans = size(data,1);

for c = 1 : n_chans
    
    channel = c;
    
    % Perform time-frequency decomposition to obtain 
    % the power at each frequency throughout time
    if strcmp(tf_method,'wavelet')
        
        [ ~ , freq_vector, ~ , power ] = ...
            tf_decomposition_wavelet(data, c, ...
            f_range, n_freq, fs_eeg, 1);       
   
    elseif strcmp(tf_method,'welch')
        
        [ freq_vector, power ] =  ...
            tf_decomposition_welch(data, channel, ...
            f_range, n_freq, fs_eeg, fs, 0);
        
    end
    
    power_all(:,:,c) = power;
        
end