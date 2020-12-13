function [features_delayed] = convolve_features(features,fs,delay_range, ...
    kern_seconds)
% 
%   function [feature_delay] = convolve_features(feature,fs,delay_range,
%   kern_seconds) convolves each feature in variable 'feature' with a set  
%   of HRF functions, characterized by a set of overshoot delays specified
%   in 'delay_range'
%
%   INPUTS:
%
%   feature         the features to be convolved  
%   fs              the sampling frequency of the features
%   delay_range     the vector specifying the delays (sec)
%   kern_seconds    the size of the hrf kernel to be used for 
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

% Read size of input data
siz = size(features);
n_pnts = siz(1);
n_features = prod(siz(2:end));
n_delays = length(delay_range);

% ------------------------------------------------------------
% Preallocate matrices 
% ------------------------------------------------------------     

% Allocate new features matrix, by adding the delays dimension 
features_delayed = zeros(n_pnts, n_features, n_delays);

% Reshape the features matrix to have only 2 dimensions 
features = reshape(features, [n_pnts n_features]);

% Allocate matrix of hrfs 
kern_samples = kern_seconds*fs + 1;
hrf = zeros(kern_samples, length(delay_range));

% Assign hrf basis function struct, in 
% the format required by spm_Volterra 
xBF.dt = 1/fs;
xBF.name = 'hrf';
xBF.length = kern_seconds;
xBF.order = 1;

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
    p = [overshoot_delay 16*s 1*s 1*s 6 0 kern_seconds];
    
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
    bf = xBF.bf;

    % Perform convolution between each 
    % feature and each HRF delay 
    parfor f = 1 : n_features 
        
        % Create eeg predictor struct, in 
        % the format required by spm_Volterra
        P = []; P.name = {'feature'};

        % Convolution in the time-domain, with the same size as
        % the original signal (first and last samples of
        % convolution are cut out to match the intended size)
        P.u = features(:, f); 
        features_delayed(:, f, n)= spm_Volterra(P, bf);    

    end
    
    n = n + 1;
    clear xBF.bf
    
end

features_delayed = reshape(features_delayed, ...
    [n_pnts siz(2:end) n_delays]);

% ------------------------------------------------------------
% Plot the HRFs   
% ------------------------------------------------------------     

% Assign time vector with fs 250 Hz 
time_hrf = 0 : 1/fs : kern_seconds;

% Plot hrfs 
figure('Name','Hemodynamic Response Functions (HRFs)')

for d = 1 : length(delay_range)
    plot(time_hrf, hrf(:,d)); hold on 
end

title('Hemodynamic Response Functions (HRFs)','FontSize',16);
xlabel('Time (s)','FontSize',16); ylabel('Amplitude','FontSize',16);
%legend("10 seconds", "8 seconds", "6 seconds",...
%    "5 seconds", "4 seconds", "2 seconds",'FontSize',14)
grid on 


end
