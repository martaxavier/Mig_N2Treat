function [feature_delay] = convolve_features(feature,fs,delay_range,size_kern)
% 
%   function [feature_delay] = convolve_features(feature,fs,delay_range,
%   size_kern) convolves each feature in variable 'feature' with a set  
%   of HRF functions, characterized by a set of overshoot delays specified
%   in 'delay_range'
%
%   INPUTS:
%
%   feature         the features to be convolved  
%                   (time x chans x freq bands) 
%   fs              the sampling frequency of the features
%   delay_range     the vector specifying the delays (sec)
%   size_kern       the size of the hrf kernel to be used for 
%                    convolution (sec)
%
%   OUTPUTS:
%
%   feature_delay   the features after being convolved
%                   (time x chans x delays x freq bands)
%

% ------------------------------------------------------------
% Read input information
% ------------------------------------------------------------     

% Extract data from pred
n_chans = size(feature,2);               % number of channels
n_pnts = size(feature,1);                % number of data points

% ------------------------------------------------------------
% Preallocate matrices 
% ------------------------------------------------------------     

% Allocate feature_delay matrix for all frequency bands
% time goes through rows, channels go through cols, 
% delays go through depth, predictor go through 4D
feature_delay = zeros(n_pnts, n_chans, ...
    length(delay_range),size(feature,3)); 

% Allocate matrix of hrfs 
n_kern = size_kern*fs + 1;
hrf = zeros(n_kern,length(delay_range));

% Assign hrf basis function struct, in 
% the format required by spm_Volterra 
xBF.dt = 1/fs;
xBF.name = 'hrf';
xBF.length = size_kern;
xBF.order = 1;

% Create eeg predictor struct, in 
% the format required by spm_Volterra
P.name = {'predictor'};

% ------------------------------------------------------------
% Convolve matrices with hrf kernel
% ------------------------------------------------------------     

n = 1;

% Go through overshoot delays 
for overshoot_delay = delay_range

    % Assign parameters of the response function.
    % for overshoot delay (p1) = 6, undershoot delay (p2) = 16,
    % dispersion of response (p3) = 1, dispersion of undershoot (p4) = 1,
    % ratio of response to undershoot (p5) = 6, onset (p6) = 0, 
    % length of kernel (p7) = 32 (default)
    % maintain a linear relation between the parameters’ values
    % of the five variants HRFs and the canonical HR
    s = overshoot_delay/6; % scale factor 
    p = [overshoot_delay 16*s 1*s 1*s 6 0 size_kern];
    
    % Assign scan repetition time
    % this should result in a function 
    % with the same temporal resolution 
    % as the original dataset (0.004 s)
    rt = 1/fs;

    % Build hrf function for current overshoot delay 
    % each column corresponds to hrf with a given overshoot delay
    hrf(:,n) = spm_hrf(rt, p);

    % Normalize overshoot 
    hrf(:,n) = hrf(:,n)./max(hrf(:,n));
    
    xBF.bf = hrf(:,n);

    % Perform convolution between predictor 
    % (for each channel) and hrf (for each delay)
    for channel = 1 : n_chans
        
        for i = 1 : size(feature,3)
            
             % Convolution in the time-domain, with the same size as
             % the original signal (first and last samples of
             % convolution are cut out to match the intended size)
             P.u = feature(:,channel,i); 
             feature_delay(:,channel,n,i)= spm_Volterra(P, xBF.bf);
             clear P.u;
            
        end       
        
    end
    
    n = n + 1;
    clear xBF.bf
    
end

% ------------------------------------------------------------
% Plot the hrfs  
% ------------------------------------------------------------     

% Assign time vector with fs 250 Hz 
time_hrf = 0 : 1/fs : size_kern;

% Plot hrfs 
figure('Name','Hemodynamic Response Functions (HRFs)')

for d = 1 : length(delay_range)
    plot(time_hrf, hrf(:,d)); hold on 
end

title('Hemodynamic Response Functions (HRFs)','FontSize',16);
xlabel('Time (s)','FontSize',16); ylabel('Amplitude','FontSize',16);
legend("10 seconds", "8 seconds", "6 seconds",...
    "5 seconds", "4 seconds", "2 seconds",'FontSize',14)
grid on 


end
