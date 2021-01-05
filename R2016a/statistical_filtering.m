function [conspec_sig, decision, p_thresh_corrected] = ...
    statistical_filtering(data, f_range, n_freq, n_wins_welch, ...
    win_seconds, fs_data, fs_cspec, metric, conspec, surrogate_method, ...
    n_surrogates)

% Define some parameters 
% for the surrogate analysis 
win_cut_seconds = 10;
p_thresh_surrogate = 0.05;

% Read size of input data 
n_chans = size(data,1);
n_pnts = size(data,2);

% Get metric parameters
get_metric_pars

% Obtain number of time-points of the output 
n_pnts_cspec = round((n_pnts - 1)*(fs_cspec/fs_data) + 1);

% Pre-allocate matrix of all connectivity surrogates 
conspec_surr = zeros(n_pnts_cspec, n_chans, ...
    n_chans, n_bands, n_surrogates);

bands = bands;

%---------------------------------------------------------    
% Generate connectivity surrogates  
%---------------------------------------------------------
parfor s = 1 : n_surrogates
    
    data_par = data; 
    bands_par = bands;
    
    % Generate surrogate time-series 
    switch surrogate_method
        
        case 'block_shuffle'
         
            % Create window in which the 
            % time-series is restricted 
            % to be cut 
            mid = floor(n_pnts/2);
            win_cut_samples = win_cut_seconds * fs_data + 1;
            win_cut = mid - floor(win_cut_samples/2) : ...
                mid + floor(win_cut_samples/2);
            
            % Choose randomly one point from that window 
            cut = round(win_cut(1) + rand(1)*...
                (win_cut(end) - win_cut(1)));
            
            % Cut the signal at that point and switch  
            % the two resulting halves 
            idxs_surr = [cut + 1 : n_pnts 1 : cut];
            data_surr = data_par(:,idxs_surr);
      
    end
    
    %-----------------------------------------------------    
    % TF decomposition to extract the cross-spectrum  
    %-----------------------------------------------------

    % Compute the cross-spectrum at each pair of signals, 
    % throughout time 
    [cross_spectrum_surr, f_vector] = tf_analysis_cross_spectrum...
        (data_surr, f_range, n_freq, n_wins_welch, win_seconds, ...
        fs_data, fs_cspec);

    %-----------------------------------------------------   
    % Compute connectivity feature 
    %-----------------------------------------------------

    % Compute connectivity metric 
    conspec_surr_par = compute_connectivity_metric...
        (cross_spectrum_surr, metric);

    % Average connectivity for each frequency band 
    conspec_surr_par = average_frequency...
        (conspec_surr_par, f_vector, bands_par);

    % Turn lower triangular into symmetric matrix 
    conspec_surr_par = tril2symmetric(conspec_surr_par);       
    conspec_surr(:, :, :, :, s) = conspec_surr_par;
    
end
        
%---------------------------------------------------------    
% Statistical filtering 
%---------------------------------------------------------

% p_value is the poportion of surrogates that 
% are higher than the original connectivity 
is_greater = conspec_surr > repmat(conspec, ...
    1, 1, 1, 1, n_surrogates);
p_values = sum(is_greater, 5) ./ n_surrogates;

% Correct for multiple comparisons, using the 
% FDR correction 
[decision, p_thresh_corrected, ~, ~] = ...
    fdr_bh(p_values,p_thresh_surrogate,'pdep','no');

% Set to zero connectivity values that 
% didn't survive the significance test
conspec_sig = zeros(size(conspec));
conspec_sig(decision) = conspec(decision); 