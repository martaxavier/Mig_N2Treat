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

[con_metrics, net_metrics] = get_con_net_metrics(metrics);

for m = 1 : length(con_metrics)
    
    % Define current metric
    con_metric = con_metrics{m};
    metric = con_metric;
    get_metric_pars;
    
    %------------------------------------------------------------    
    % Go through subjects
    %------------------------------------------------------------
    for s = 1 : length(subjects)

        % Define current subject 
        subject = subjects(s);
        
        disp(char(strcat('Computing connectivity features for', ... 
            {' '}, subject,', metric',{' '}, metric, ' ...')));

        %---------------------------------------------------------    
        % Load data 
        %--------------------------------------------------------- 

        % Load eeg dataset of specified subject 
        load(char(fullfile(path_data_in(s),data_in)));
        
        % Save chanlocs structure if non existent 
        if ~exist(char(fullfile('PARS','chanlocs.mat')),'file')
            save(char(fullfile('PARS','chanlocs.mat')),'chanlocs');
        end

        first_last = dlmread(char(fullfile(path_markers_in(s),markers_in)));
        
        my_file = char(fullfile(path_markers_in(s), markers_sub_task_in));
        if exist(my_file,'file')
            sub_task_design = dlmread(char(my_file));
        end
        
        % Assign first and last 
        % EEG samples
        first_eeg = first_last(1); 
        last_eeg = first_last(end);
        
        % Create output directories if non-existent 
        if ~exist(char(path_data_out(s)), 'dir'); mkdir(char(path_data_out(s))); end
        if ~exist(char(path_img_out(s)), 'dir'); mkdir(char(path_img_out(s))); end    

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
            (data, [f_min f_max], n_freq, n_wins_welch, ...
            tf_sliding_window_seconds, fs_eeg, fs);
        
        % Update number of time-points 
         n_pnts = size(cross_spectrum,2);
        
        %---------------------------------------------------------    
        % Compute connectivity feature 
        %---------------------------------------------------------
        
        % Compute connectivity metric 
        conspec = compute_connectivity_metric(cross_spectrum, con_metric);
        
        % Average connectivity for each frequency band 
        conspec_avg = average_frequency(conspec, f_vector, bands);
        
         %---------------------------------------------------------    
         % Plots and report of connectivity spectrum 
         %--------------------------------------------------------- 

%         if flag.report ~= 0
%             report_connectivity_spectrum;
%         end
    
        %---------------------------------------------------------    
        % Statistical filtering  
        %---------------------------------------------------------

        % Compute surrogates of the connectivity matrix to perform 
        % statistical filtering and remove non-significant connections 
        [conspec_sig, decision_sig, p_thresh_corrected] = ...
            statistical_filtering(data, [f_min f_max], n_freq, ...
            n_wins_welch, tf_sliding_window_seconds, fs_eeg, fs, ...
            con_metric, conspec_avg, surrogate_method, n_surrogates); 
        
        %---------------------------------------------------------    
        % Topological filtering  
        %---------------------------------------------------------

        % Perform topological filtering of the connectivity matrix
        % using Orthogonalized Minimal Spanning Trees 
        [conspec_topo, decision_topo] = topological_filtering(conspec_sig, ...
            decision); 

        % Save decision matrices and filtered connectivity matrix
        % in the output directory
        decision_sig_out = strcat('decision_sig_',con_metric,'.mat');
        decision_topo_out = strcat('decision_topo_',con_metric,'.mat');
        conspec_out = strcat('conspec_filtered_',con_metric,'.mat');
        save(char(fullfile(path_data_out(s), decision_sig_out)), 'decision_sig');
        save(char(fullfile(path_data_out(s), decision_topo_out)), 'decision_topo');
        save(char(fullfile(path_data_out(s), conspec_out)), 'conspec_topo');
                    
        % Save p-value of connectivity significance, corrected by FDR
        p_out = strcat('p_thresh_corrected_', con_metric, '.mat');
            save(char(fullfile(path_data_out(s), p_out)), 'p_thresh_corrected');

        %---------------------------------------------------------    
        % Compute network measure 
        %---------------------------------------------------------
        
        for n = 1 : length(net_metrics{:,m})
            
            net_metric = net_metrics{m};
            net_metric = net_metric(n);
            full_metric = strcat(con_metric,'_',net_metric);
            metric = full_metric;
            get_metric_pars;
            
            eeg_features = compute_network_metric(conspec_topo, ...
                net_metric); 
        
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
            % Plot network features 
            %--------------------------------------------------------- 
            
            if flag.report ~= 0     
                report_network_features;
            end
            
            %---------------------------------------------------------    
            % Write feature files 
            %--------------------------------------------------------- 

            % Build final feature matrices to be written in file        
            % Build final feature matrices to be written in file        
            eeg_features_norm = reshape(eeg_features_norm, ...
                [n_pnts, numel(eeg_features_norm(1, :, :, :))]);

            % Delete the channel pairs that correspond to redundant 
            % channel pairs 
            if length(squeeze(dim)) == 4
                eeg_features_norm = eeg_features_norm(:, find(repmat...
                    (tril(ones(n_chans), -1), [1 1 n_bands])));
            end

            if ~isempty(eeg_shift)

                eeg_features_delayed_norm = reshape(...
                    eeg_features_delayed_norm, [n_pnts, ...
                    numel(eeg_features_delayed_norm(1, :, :, :, :))]);

            end

            % Specify file names      
            feature_out = char(strcat(eeg_metric, '_', data_out));

            switch eeg_shift

                case 'conv'

                    feature_out_delayed = char(strcat(eeg_metric,'_', ...
                        data_out_conv));

                case 'delay'

                    feature_out_delayed = char(strcat(eeg_metric,'_', ...
                        data_out_delay));

            end 

            % Write EEG feature files
            dlmwrite(char(fullfile(path_data_out(s), feature_out)), ...
                eeg_features_norm);

            if ~isempty(eeg_shift)

                dlmwrite(char(fullfile(path_data_out(s),feature_out_delayed)), ...
                    eeg_features_delayed_norm);

            end

        end % finish looping through network metrics 
        
    end % finish looping through subjects 

end % finish looping through connectivity metrics 

%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% function get_con_net_metrics(metrics)          
% ============================================================

function [all_con_metrics, all_net_metrics] = ...
    get_con_net_metrics(metrics)

% Pre-allocate strings of that list the 
% connectivity and network metrics 
con_metrics = strings(length(metrics), 1);
net_metrics = strings(length(metrics), 1);

for m = 1 : length(metrics)
    con_metrics(m) = extractBefore(metrics(m), '_');
    net_metrics(m) = extractAfter(metrics(m), '_');
end

[con_metrics, ~, c] = unique(con_metrics);

% Pre-allocate 
all_con_metrics = cell(length(unique(c)), 1);
all_net_metrics = cell(1, length(unique(c)));

for m = 1 : length(unique(c))
    all_con_metrics{m} = con_metrics(m);
    all_net_metrics{:,m} = net_metrics(c == m);
end

end