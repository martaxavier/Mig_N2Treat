% Computes the EEG connectivity features to be used as the model's
% design matrix X
%
% 	The method here implemented is as follows: 
%
%       1)  Time-frequency decomposition of the EEG data at each
%           channel to obtain the cross-spectrum at each pair of 
%           channels - and simultaneous downsampling the the EEG
%           data to the analysis's sampling frequency (welch method)
%
%       2)  Computation of the specified connectivity metrics from
%           the cross-spectra and averaging across frequency bands
%
%       3)  Convolution of all features with a range of HRFs/ 
%           shifting of all the features with a range of delays 
%
%       4)  Prunning of all features according to match the 
%           beginning and end of the simultaneous BOLD acquisition/
%           according to the current task design
%
%       6)  Normalization of all EEG features to have 0 mean and 
%           standard deviation 1
%

% Intermediate fs (Hz)
fs = fs_analysis;          

%------------------------------------------------------------    
% Go through metrics
%------------------------------------------------------------
for m = 1 : length(metrics)
    
    % Define current metric
    metric = metrics(m);
    
    % Get parameters of the 
    % current metric
    get_metric_pars
    
    %------------------------------------------------------------    
    % Go through subjects
    %------------------------------------------------------------
    for s = 1 : length(subjects)

        % Define current subject 
        subject = subjects(s);
        
        disp(strcat('Computing connectivity features for', ... 
            " ", subject,', metric'," ", metric, ' ...'));

        %---------------------------------------------------------    
        % Load data 
        %--------------------------------------------------------- 

        % Load eeg dataset of specified subject 
        load(fullfile(path_data_in(s),data_in));
        
        % Save chanlocs structure if non existent 
        if ~exist(fullfile('PARS','chanlocs.mat'),'file')
            save(fullfile('PARS','chanlocs.mat'),'chanlocs');
        end

        first_last = dlmread(fullfile(path_markers_in(s),markers_in));
        
        my_file = fullfile(path_markers_in(s), markers_sub_task_in);
        if exist(my_file,'file')
            sub_task_design = dlmread(my_file);
        end
        
        % Assign first and last 
        % EEG samples
        first_eeg = first_last(1); 
        last_eeg = first_last(end);
        
        % Create output directories if non-existent 
        if ~exist(path_data_out(s), 'dir'); mkdir(path_data_out(s)); end
        if ~exist(path_img_out(s), 'dir'); mkdir(path_img_out(s)); end    

        % Extract EEG data from dataset
        data = EEG.data;                 

        % First and last samples in the new fs
        first = ceil(first_eeg*fs/fs_eeg);
        last = ceil(last_eeg*fs/fs_eeg); 

        %---------------------------------------------------------    
        % TF decomposition to extract the cross-spectrum  
        %---------------------------------------------------------
        
        % Compute the cross-spectrum at each pair of signals, 
        % throughout time 
        [cross_spectrum, f_vector] = tf_analysis_cross_spectrum...
            (data, [f_min f_max], n_freq, tf_sliding_window_seconds, ...
            fs_eeg, fs);
        
        % Update number of time-points 
         n_pnts = size(cross_spectrum,2);
        
        %---------------------------------------------------------    
        % Compute connectivity feature 
        %---------------------------------------------------------
        
        % Estimate number of samples of the TF sliding window 
         tf_sliding_window_samples = tf_sliding_window_seconds...
             *fs_eeg + 1;
                        
        % Compute connectivity matric 
        [conspec, pvalues] = compute_connectivity_metric...
            (cross_spectrum, metric, tf_sliding_window_samples);
        
        % Average connectivity for each frequency band 
        eeg_features = average_frequency(conspec, f_vector, bands);
        
        %---------------------------------------------------------    
        % Remove outliers 
        %---------------------------------------------------------

        % Remove feature time-points that are 10 times greater
        % than the features' standard deviation 
        
        % Reshape the feature matrix to have 2 dimensions
        siz = size(eeg_features);
        eeg_features = reshape(eeg_features,[n_pnts, ...
            numel(eeg_features(1, :, :, :))]);
        
        % Define the outlier threshold for each feature 
        out_thresh = repmat(10*std(eeg_features),[n_pnts 1]);
        eeg_features(eeg_features>out_thresh) = ...
            out_thresh(eeg_features>out_thresh);
        
        % Reshape feature matrix into its original dimension 
        eeg_features = reshape(eeg_features, siz);
        
        %---------------------------------------------------------    
        % Mirror padd features before convolution  
        %---------------------------------------------------------        
        
        % Mirror padd the features at a length equal
        % to the convolution kernal size + 1 (seconds)
        eeg_features = eeg_features(first:last, :, :, :);

        % Padd features in the pre-direction 
        padsize = max(first - 1, hrf_kernel_seconds*fs);
        eeg_features = padarray(eeg_features, ...
            padsize, 'symmetric', 'pre');

        % Padd features in the post-direction
        padsize = max(n_pnts - last, hrf_kernel_seconds*fs);
        eeg_features = padarray(eeg_features, ...
            padsize, 'symmetric', 'post');
        
        %---------------------------------------------------------    
        % Convolve/delay features 
        %--------------------------------------------------------- 
        
        switch eeg_shift
            
            case 'conv'
                
                % Specify label for plotting 
                plotting_shift = 'convolved';
                
                % Convolve features with the specified family of HRFs 
                eeg_features_delayed = convolve_features(eeg_features, ...
                    fs, delays, hrf_kernel_seconds);
                
                % Permute resulting matrix to have bands at the end
                siz = size(eeg_features_delayed);
                eeg_features_delayed = permute(eeg_features_delayed, ...
                    [(1:length(siz(1:end-2))) length(siz) ...
                    length(siz(1 : end-1))]);
                
                % Prune the features to match the BOLD acquisition
                % times
                eeg_features_delayed = eeg_features_delayed...
                    (first:last, :, :, :, :);                
                
            case 'delay'
                
                % Specify label for plotting 
                plotting_shift = 'delayed';
                
                % Shift features by the specified range of delays 
                eeg_features_delayed = delay_features(eeg_features, ...
                    fs, delays, [first last]);
                
                % Prune the features to match the BOLD acquisition
                % times
                eeg_features_delayed = eeg_features_delayed...
                    (first:last, :, :, :, :);
             
        end
         
        % Prune the original features to match the BOLD 
        % acquisition times 
        eeg_features = eeg_features(first:last, :, :, :, :);
        
        % Update number of time-points 
        n_pnts = size(eeg_features,1);
        time = 0 : 1/fs : (n_pnts-1)/fs;
            
        %---------------------------------------------------------    
        % Prune according to the task design (case sub_task)
        %---------------------------------------------------------
        
        if ~isempty(sub_task)
            
            % Prune the sub_task_design to match the duration 
            % of the task 
            sub_task_design = sub_task_design(first : last);
            
            % Prune the original features according 
            % to the sub-task design 
            eeg_features = eeg_features(logical...
                (sub_task_design), :, :, :);
            
            if ~isempty(eeg_shift)
                
                % Prune the delayed features according to the 
                % sub-task design 
                eeg_features_delayed = eeg_features_delayed...
                    (logical(sub_task_design), :, :, :, :);
                
            end
            
            % Update number of time-points 
            n_pnts = size(eeg_features,1);
            time = 0 : 1/fs : (n_pnts-1)/fs;
            
        end
        
        %---------------------------------------------------------    
        % Normalize features (0 mean, 1 std) 
        %--------------------------------------------------------- 
        
        % Normalize (subtract the mean, divide by the std)
        eeg_features_norm = zscore(eeg_features);
        
        if ~isempty(eeg_shift)
            
            % Normalize by first reshaping features into a 
            % 2 dimensional array
            eeg_features_delayed_norm  = zscore(reshape(...
                eeg_features_delayed, [n_pnts numel(...
                eeg_features_delayed(1, :, :, :, :))]));
            
            % Reshape features into its original dimensions
            eeg_features_delayed_norm = reshape(...
                eeg_features_delayed_norm, size(...
                eeg_features_delayed));
            
        end
        
        %---------------------------------------------------------    
        % Plots and report (before/after convolution)
        %--------------------------------------------------------- 
        
        if flag.report ~= 0
            report_connectivity_features;
        end
        
        %---------------------------------------------------------    
        % Write feature files 
        %--------------------------------------------------------- 
        
        % Build final feature matrices to be written in file        
        eeg_features_norm = reshape(eeg_features_norm, ...
            [n_pnts, numel(eeg_features_norm(1, :, :, :))]);
        
        % Delete the channel pairs that correspond to redundant 
        % channel pairs 
        if size(squeeze(dim)) == 4
            eeg_features_norm = eeg_features_norm(:, find(repmat...
                (tril(ones(n_chans), -1), [1 1 n_bands])));
        end
        
        if ~isempty(eeg_shift)
            
            eeg_features_delayed_norm = reshape(...
                eeg_features_delayed_norm, [n_pnts, ...
                numel(eeg_features_delayed_norm(1, :, :, :, :))]);
            
        end
        
        % Specify file names      
        feature_out = strcat(eeg_metric, '_', data_out);

        switch eeg_shift

            case 'conv'

                feature_out_delayed = strcat(eeg_metric,'_', ...
                    data_out_conv);

            case 'delay'

                feature_out_delayed = strcat(eeg_metric,'_', ...
                    data_out_delay);
                
        end 

        % Write EEG feature files
        dlmwrite(fullfile(path_data_out(s), feature_out), ...
            eeg_features_norm);
        
        if ~isempty(eeg_shift)
            
            dlmwrite(fullfile(path_data_out(s),feature_out_delayed), ...
                eeg_features_delayed_norm);
            
        end
        
    end % finish looping through subjects 
    
end % finish looping through metrics 
