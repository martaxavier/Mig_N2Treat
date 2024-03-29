% Computes the EEG power features to be used as the model's design
% matrix X
%
% 	The method here implemented is as follows: 
%
%       1)  Time-frequency decomposition of the EEG data at each
%           channel
%
%       2)  Computation of the specified EEG features from the 
%           time-frequency spectrum of the EEG data 
%
%       3)  Convolution of all features with a range of HRFs/
%           shifting of all features with a Frange of delays 
%
%       4)  Psening of all features according to the beginning 
%           and end of the simultaneous BOLD acquisition/ according
%           to the current task design 
%
%       5)  Downsampling of all EEG features to the analysis's fs
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
    
    %------------------------------------------------------------    
    % Go through subjects
    %------------------------------------------------------------
    for s = 1 : length(subjects)

        % Define current subject 
        subject = subjects(s);
        
        disp(strcat('Computing power features for', ... 
            " ", subject,', metric'," ", metric, ' ...'));

        %---------------------------------------------------------    
        % Load data 
        %--------------------------------------------------------- 

        % Create output directories if non-existent 
        if ~exist(path_data_out(s, se), 'dir'); mkdir(path_data_out(s, se)); end
        if ~exist(path_img_out(s, se), 'dir'); mkdir(path_img_out(s, se)); end   
        
        % Load eeg dataset of specified subject 
        load(fullfile(path_data_in(s, se), data_in));

        first_last = dlmread(fullfile(path_markers_in(s, se), markers_in));
        
        my_file = fullfile(path_markers_in(s, se), markers_sub_task_in);
        if exist(my_file, 'file')
            sub_task_design = dlmread(my_file);
        end
        
        % Assign first and last 
        % EEG samples
        first_eeg = first_last(1); 
        last_eeg = first_last(end); 

        % Extract EEG data from dataset
        data = EEG.data;                 

        % First and last samples in the new fs
        first = ceil(first_eeg*fs/fs_eeg);
        last = ceil(last_eeg*fs/fs_eeg); 

        %---------------------------------------------------------    
        % TF decomposition to extract power 
        %---------------------------------------------------------
        
        disp('... performing time-frequency decomposition ...');
                    
        [power, f_vector] = tf_analysis_power_spectrum...
            (data, [f_min f_max], n_freq, tf_method, ...
            tf_wavelet_kernel_seconds, tf_sliding_win_seconds, ...
            n_wins_welch, fs_eeg, fs);   
 
        n_pnts = size(power, 1);

        %---------------------------------------------------------    
        % Compute power metric
        %---------------------------------------------------------
        
        get_metric_pars

        % Lin comb of power 
        % across frequency bands 
        if n_bands > 1 

            % Average power across frequency bands
            eeg_features = average_frequency(power, ...
                f_vector, bands, 3);

        % Total power 
        elseif contains(metric, 'tp')

            % Compute total power 
            eeg_features = squeeze(sum(power, 3));

        % Root mean square frequency
        elseif contains(metric, 'rmsf')

            % Compute root mean square frequency
            eeg_features = squeeze(sqrt(sum(permute(repmat...
                ((f_vector' .^ 2), 1, n_pnts, n_chans), ...
                [2 3 1]) .* power, 3)));
            
        end

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
        out_thresh = repmat(10*std(eeg_features), [n_pnts 1]);
        eeg_features(eeg_features > out_thresh) = ...
            out_thresh(eeg_features > out_thresh);
        
        % Reshape feature matrix into its original dimension 
        eeg_features = reshape(eeg_features, siz);

        %---------------------------------------------------------    
        % Assign current first and last indices 
        %---------------------------------------------------------

        % Assign current first and last indices according
        % to the current sampling frequency of the data 
        if strcmp(tf_method,'welch')
            
            first_current = first; 
            last_current = last;
            fs_current = fs;
            
        elseif strcmp(tf_method,'wavelet')
            
            first_current = first_eeg; 
            last_current = last_eeg;
            fs_current = fs_eeg;
            
        end

        %---------------------------------------------------------    
        % Convolution with HRF or delay
        %---------------------------------------------------------

        disp('... convolving with a range of HRFs ...');
        
        switch eeg_shift 

            case 'conv'

                % Specify label for plotting 
                plotting_shift = 'convolved';
                
                % Convolve features with the specified family of HRFs 
                eeg_features_delayed = convolve_features_fast(eeg_features, ...
                    fs_current, delays, hrf_kernel_seconds);

                % Permute resulting matrix to have bands at the end
                if n_bands > 1
                    perm = 1 : ndims(squeeze(eeg_features_delayed));
                    eeg_features_delayed = permute...
                        (eeg_features_delayed, ...
                        [perm(1:end-2) perm(end) perm(end-1)]);
                end
                
                % Psee the features to match the BOLD acquisition
                % times
                eeg_features_delayed = eeg_features_delayed...
                    (first_current:last_current, :, :, :);

            case 'delay'

                % Specify label for plotting 
                plotting_shift = 'delayed';
                
                % Shift features by the specified range of delays 
                eeg_features_delayed = delay_features(eeg_features, ...
                    fs_current, delays, [first_current,last_current]); 
                
                % Psee the features to match the BOLD acquisition
                % times
                eeg_features_delayed = eeg_features_delayed...
                   (first_current:last_current, :, :, :, :);

        end
        
        % Prune the original features to match the BOLD acquisition times
        eeg_features = eeg_features(first_current : last_current, :, :);

        %---------------------------------------------------------    
        % Prune according to the task design (case sub_task)
        %---------------------------------------------------------
        
        if ~isempty(sub_task)
            
            % Psee the sub_task_design to match the duration 
            % of the task 
            sub_task_design = sub_task_design(first_current : ...
                last_current);
            
            % Psee the original features according to the
            % sub-task design 
            eeg_features = eeg_features(logical(sub_task_design), ...
                :, :);
            
            if ~isempty(eeg_shift)
                
                % Psee the delayed features according to the 
                % sub-task design 
                eeg_features_delayed = eeg_features_delayed...
                    (logical(sub_task_design), :, :, :);
                
            end
            
        end

        %---------------------------------------------------------    
        % Downsample to intermediate frequency
        %---------------------------------------------------------

        disp('... downsampling to intermeadiate frequency ...');
        
        % Keep a version of the features in  
        % the original EEG sampling frequency 
        eeg_features_eeg_fs = eeg_features;
        
        % Update the number of time-points and 
        % and the time vector according to the 
        % current sampling frequency
        n_pnts_current = size(eeg_features,1);
        time_current = 0 : 1/fs_current : ...
            (n_pnts_current - 1)/fs_current;
        
        switch tf_method

            case 'welch'

               % If the Welch method was used, the 
               % data was already downsampled 
               n_pnts = n_pnts_current;
               time = time_current;

            case 'wavelet'
                
                % If Morlet wavelet decomposition was used, the 
                % data is yet in the original sampling frequency 
                time = 0 : 1/fs : (n_pnts_current-1)/fs_current; 
                n_pnts = length(time);

                % Use spline interpolation to downsample features 
                eeg_features = permute(eeg_features, [3 2 1]);         
                eeg_features = spline(time_current, eeg_features, time);         
                eeg_features = permute(eeg_features, [3 2 1]);                

                if ~isempty(eeg_shift)
                    
                    eeg_features_delayed = ...
                        permute(eeg_features_delayed, [4 3 2 1]);    
                    eeg_features_delayed = ...
                        spline(time_current, eeg_features_delayed, time);
                    eeg_features_delayed = ...
                        permute(eeg_features_delayed, [4 3 2 1]);  
                    
                end
                
        end 

        %---------------------------------------------------------    
        % Normalize features (0 mean, 1 std)
        %---------------------------------------------------------

        disp('... normalizing features ...');
        
        % Normalize every feature so as to have
        % zero mean and standard deviation one
        eeg_features_norm = zscore(eeg_features);

        if ~isempty(eeg_shift)
       
            eeg_features_delayed_norm = ...
                zscore(reshape(eeg_features_delayed, [n_pnts, ...
                numel(eeg_features_delayed(1, :, :, :))]));
            eeg_features_delayed_norm = ...
                reshape(eeg_features_delayed_norm, ...
                size(eeg_features_delayed));
            
        end

        %---------------------------------------------------------    
        % Plot features (after normalization)
        %---------------------------------------------------------
        
        if ~isempty(eeg_shift) && flag.report ~= 0
            
            disp('... generating and saving plots of features ...');
                    
            for c = 1 : n_chans

                for b = 1 : n_bands

                    figure_name = strcat(id_bands(b), ...
                        ' EEG features of'," ",subject, ...
                        ', channel'," ", num2str(c,'%01d')); 
                    figure('Name',figure_name)

                    plot(time, eeg_features_norm(:,c,b)); hold on; 
                    plot(time, eeg_features_delayed_norm(:,c, ...
                        plotting_delay,b));

                    title(strcat(id_bands(b), ' feature of'," ", ...
                        subject,'channel'," ", num2str(c,'%01d'), ...
                        ' (demeaned, normalized)')); 
                    xlabel('Time(s)'); ylabel('Power (\muV^2)');
                    legend(strcat(id_bands(b),' feature'), ...
                        strcat(id_bands(b),' feature,', ...
                        plotting_shift,', delay'," " , ...
                        num2str(delays(plotting_delay)), 's'));

                    img_out = strcat(id_bands(b), num2str(delays ...
                        (plotting_delay)), 'sec', num2str(c), 'chan.png');
                    saveas(gcf,fullfile(path_img_out(s, se), img_out));

                end % looping through bands

            end % looping through channels
            
        end
        
        %---------------------------------------------------------    
        % Write feature files 
        %---------------------------------------------------------

        % Build final feature matrix to be written in file
        eeg_features_norm = reshape(eeg_features_norm, ...
            [n_pnts,numel(eeg_features_norm(1,:,:))]);
        eeg_features_eeg_fs = reshape(eeg_features_eeg_fs, ...
            [size(eeg_features_eeg_fs,1), ...
            numel(eeg_features_eeg_fs(1,:,:))]);
        
        if ~isempty(eeg_shift)
            
            eeg_features_delayed_norm = ...
                reshape(eeg_features_delayed_norm, ...
                [n_pnts,numel(eeg_features_delayed_norm(1,:,:,:))]);
            
        end
        
        % Specify file names
        switch eeg_shift

            case 'conv'

                feature_out = strcat(eeg_metric,'_',data_out);
                feature_out_eeg_fs = ...
                    strcat(eeg_metric,'_',data_out_eeg_fs);
                feature_out_delay = ...
                    strcat(eeg_metric,'_',data_out_conv);

            case 'delay'

                feature_out = strcat(eeg_metric,'_',data_out);
                feature_out_eeg_fs = ...
                    strcat(eeg_metric,'_',data_out_eeg_fs);
                feature_out_delay = ...
                    strcat(eeg_metric,'_',data_out_delay);
                
            case ''
                
                feature_out = strcat(eeg_metric,'_',data_out);
                feature_out_eeg_fs = ...
                    strcat(eeg_metric,'_',data_out_eeg_fs);
                
        end 

        % Write EEG feature files
        dlmwrite(fullfile(path_data_out(s, se), ...
            feature_out), eeg_features_norm);

        %dlmwrite(fullfile(path_data_out(s), ...
        %    feature_out_eeg_fs),eeg_features_eeg_fs);
        
        if ~isempty(eeg_shift)
            dlmwrite(fullfile(path_data_out(s), ...
                feature_out_delay),eeg_features_delayed_norm);
        end

    end % finish looping through subjects 

end % finish looping through metrics 