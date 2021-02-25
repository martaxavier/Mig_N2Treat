% Computes the EEG connectivity features to be used as the model's
% design matrix X
%
% 	The method here implemented is as follows: 
%
%   1)  Estimation of the dynamic functional connectivity through a 
%       windowing procedure, and simultaneous downsampling of the EEG
%       data. In each window:
%           1.1)    Time-frequency decomposition of the EEG data at each 
%                   channel to obtain the expected cross-spectrum at each
%                   pair of channels, using the Welch's method
%           1.2)    Estimation of the specified connectivity metric 
%
%   2)  Averaging of the connectomes at each frequency band 
% 
%   3)  Estimation of the specified graph measure from the connectomes 
% 
%   4)  Convolution of the resulting EEG features with a range of HRFs/
%       shifting of all the features with a range of delays 
% 
%   5)  Prunning of all the features to match the beginning and end of the 
%       simultaneous BOLD acquisition/according to the current task design
% 
%      6)   Normalization of all EEG features to have 0 mean and standard 
%           deviation 1  
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
        
        disp(strcat('Computing connectivity features for', ... 
            " ", subject, ', metric', " ", metric, ' ...'));

        %---------------------------------------------------------    
        % Load data 
        %--------------------------------------------------------- 

        % Load eeg dataset of specified subject 
        load(fullfile(path_data_in(s), data_in));
        
        % Save chanlocs structure if non existent 
        if ~exist(fullfile(dataset, 'PARS', 'chanlocs.mat'), 'file')
            save(fullfile(dataset, 'PARS',' chanlocs.mat'), 'chanlocs');
        end

        first_last = dlmread(fullfile(path_markers_in(s), markers_in));
        
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
        % TF decomposition to extract the connectivity spectrum   
        %---------------------------------------------------------
        
        disp('... computing the connectomes ...'); 
        
        % Compute the cross-spectrum at each pair of signals, throughout
        % time and update vector of frequencies 
        [conspec, p_values, f_vector] = tf_analysis_con_spectrum(data, ...
            [f_min f_max], n_freq, tf_sliding_win_seconds, fs_eeg, ...
            fs, con_metric, stat_filt_method);
        
        % Update number of time-points 
         n_pnts = size(conspec,1);        
           
        %---------------------------------------------------------    
        % Statistical filtering + frequency average 
        %---------------------------------------------------------

        disp('... filtering the connectomes for significance ...');
        
        % Compute surrogates of the connectivity matrix to perform 
        % statistical filtering and remove non-significant connections 
        [conspec_sig, decision_sig, p_values, p_thresh_corrected] = ...
            statistical_filtering(data, [f_min f_max], n_freq, ...
            f_vector, tf_sliding_win_seconds, fs_eeg, fs, ...
            con_metric, conspec, p_values, stat_filt_method, ...
            surr_method, n_surrs, path_img_out(s));
        
        % Average connectivity for each frequency band in case it hasn't
        % been done yet 
        if ~strcmp(stat_filt_method, 'surrogate')
            conspec = average_frequency(conspec, f_vector, bands, 3);
            p_values = average_frequency(p_values, f_vector, bands, 3); 
            conspec_sig = average_frequency(conspec_sig, f_vector, bands, 3);
        end
  
        % Convert connectivity matrices into symmetrical matrices
        % prior to topological filtering and plotting 
        conspec = vec2tri(conspec, 2, n_chans, 'upper');
        conspec = tri2symmetric(conspec, 2, 3);
        conspec(isnan(conspec)) = 0;  
        
        % Convert p-value matrices into symmetrical matrices
        % prior to topological filtering and plotting      
        p_values = vec2tri(p_values, 2, n_chans, 'upper');
        p_values = tri2symmetric(p_values, 2, 3);
        p_values(isnan(p_values)) = 0;
        
        % Convert connectivity matrices into symmetrical matrices
        % prior to topological filtering and plotting 
        conspec_sig = vec2tri(conspec_sig, 2, n_chans, 'upper');
        conspec_sig = tri2symmetric(conspec_sig, 2, 3);
        conspec_sig(isnan(conspec_sig)) = 0;           
        
        %---------------------------------------------------------    
        % Topological filtering  
        %---------------------------------------------------------

        % Perform topological filtering of the connectivity matrix
        % using Orthogonalized Minimal Spanning Trees 
        % Only if the method for statistical filtering is not 
        % surrogate analysis 
        if strcmp(stat_filt_method, 'analytical')
            disp('... topological filtering the connectomes ...');
            [conspec_topo, decision_topo] = ...
                topological_filtering(conspec_sig); 
        else
            conspec_topo = conspec_sig; 
            decision_topo = decision_sig;
        end
 
        %---------------------------------------------------------    
        % Plots and report of the connectivity matrices   
        %--------------------------------------------------------- 
        
        % Generate images of the estimated 
        % connectivity matrices and save them 
        if flag.report ~= 0  
            disp('... generating and saving plots of connectomes ...');
            report_connectomes
        end
        
        %---------------------------------------------------------    
        % Save the connectivity matrices 
        %--------------------------------------------------------- 
        
        % Save decision matrices and filtered connectivity matrix
        % in the output directory
%         decision_sig_out = strcat('decision_sig_',con_metric,'.mat');
%         decision_topo_out = strcat('decision_topo_',con_metric,'.mat');
%         conspec_sig_out = strcat('conspec_sig_filtered_',con_metric,'.mat');
%         conspec_topo_out = strcat('conspec_topo_filtered_',con_metric,'.mat');                
%         save(fullfile(path_data_out(s), decision_sig_out), 'decision_sig');
%         save(fullfile(path_data_out(s), decision_topo_out), 'decision_topo');
%         save(fullfile(path_data_out(s), conspec_sig_out), 'conspec_sig');
%         save(fullfile(path_data_out(s), conspec_topo_out), 'conspec_topo');
                    
        % Save p-value of connectivity significance, corrected by FDR
        p_out = strcat('p_thresh_corrected_', con_metric, '.mat');
            save(fullfile(path_data_out(s), p_out), 'p_thresh_corrected');

        %---------------------------------------------------------    
        % Compute network measure 
        %---------------------------------------------------------
        
        for n = 1 : length(net_metrics{:,m})
            
            net_metric = net_metrics{m};
            net_metric = net_metric(n);
            full_metric = strcat(con_metric,'_',net_metric);
            metric = full_metric;
            get_metric_pars;
            
            disp(strcat('... computing the', " ", net_metric, ' ...'));
        
            eeg_features = compute_network_metric(conspec_topo, ...
                net_metric); 

            %---------------------------------------------------------    
            % Convolve/delay features 
            %--------------------------------------------------------- 

            disp('... convolving with a range of HRFs ...');
            
            eeg_features = eeg_features(first:last, :, :, :);
                        
            switch eeg_shift

                case 'conv'

                    % Specify label for plotting 
                    plotting_shift = 'convolved';

                    % Convolve features with the specified family of HRFs 
                    eeg_features_delayed = convolve_features_fast(eeg_features, ...
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
                disp('... generating and saving plots of features ...');
                report_network_features;
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
                eeg_features_norm = eeg_features_norm(:, logical(repmat...
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
