function [feature_delay] = delay_features(feature,fs,delay_range,samples)
% 
%   function [feature_delay] = delay_features(feature,fs,delay_range,
%   samples) shifts the input features, by a range of delays specified 
%   in 'delay_range'
%
%   INPUTS:
%
%   feature         the features to be delayed 
%                   (time x chans x freq bands) 
%   fs              the sampling frequency of the features
%   delay_range     the vector specifying the delays (sec)
%
%   OUTPUTS:
%
%   feature_delay   the features after being shifted by the delays
%                   (time x chans x delays x freq bands)
%

% ------------------------------------------------------------
% Read input information
% ------------------------------------------------------------ 

% Extract data from pred
n_chans = size(feature,2);               % number of channels

first = samples(1);                  % first eeg sample trial 
last = samples(2);                   % last eeg sample of trial 

% ------------------------------------------------------------
% Preallocate matrices 
% ------------------------------------------------------------ 

% Allocate pred_delay matrix for all frequency bands
% time goes through rows, channels go through cols, 
% delays go through depth, predictor go through 4D
feature_delay = zeros(length(first:last), n_chans,...
    length(delay_range),size(feature,3)); 

% ------------------------------------------------------------
% Delay features
% ------------------------------------------------------------ 

n = 1; 

% Go through delays 
for delay = delay_range
    
    % Update first and last sample
    % for the current delay (sec)
    current_first = first-delay*fs;
    current_last = last-delay*fs;
    
    % Go through channels 
    for channel = 1 : n_chans
        
        % Go through freq bands 
        for i = 1 : size(feature,3)
            
             feature_delay(:,channel,n,i)= ...
                 feature(current_first:current_last,...
                 channel,i);
            
        end       
        
    end
    
    n = n + 1;

end 