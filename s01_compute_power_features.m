% Computes the EEG features to be used as the model's design matrix X
%
% 	The method here implemented is as follows: 
%
%       1)  Time-frequency decomposition of the
%           EEG data at each channel
%
%       2)  Computation of the specified EEG features  
%           from the time-frequency spectrum of the 
%           EEG data 
%
%       3)  Convolution of all features with a range 
%           of HRFs/ shifting of all features with a 
%           range of delays 
%
%       4)  Prunning of all features according to the 
%           beginning and end of the simultaneous BOLD
%           acquisitio/ according to the current task 
%           design
%
%       5)  Downsampling of all EEG features to the
%           analysis' sampling frequency 
%
%       6)  Normalization of all EEG features to have 
%           0 mean and standard deviation 1
%

%------------------------------------------------------------    
% Analysis parameters
%------------------------------------------------------------

kernel_size = 32;               % size of hrf kernel (seconds)
fs = fs_analysis;               % intermediate samping frequency (Hz)

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

        % Load eeg dataset of specified subject 
        load(fullfile(path_data_in(s),data_in));
        
        % Save chanlocs structure if non existent 
        if ~exist(fullfile('PARS','chanlocs.mat'),'file')
            save(fullfile('PARS','chanlocs.mat'),'chanlocs');
        end

        first_last = dlmread(fullfile(path_markers_in(s),markers_in));
        sub_task_design = dlmread(fullfile(path_markers_in(s), ...
            markers_sub_task_in));
        
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
        % TF decomposition to extract power 
        %---------------------------------------------------------
        
        [power,f_vector] = extract_power...
            (data,[fmin fmax], n_freq, tf_method, ...
            fs_eeg, fs); 
        n_pnts = size(power,2);

        %---------------------------------------------------------    
        % Compute power metric
        %---------------------------------------------------------

        get_metric_pars
        eeg_features = zeros(n_pnts,n_chans,n_bands);

        % Lin comb of power 
        % across frequency bands 
        if n_bands > 1 

            for b = 1 : n_bands

                % Assign frequency bands indices 
                band_idx = find(f_vector>=bands(1,b) ...
                    & f_vector<bands(2,b));

                % Compute average power across current frequency band 
                eeg_features(:,:,b) = mean(power(band_idx,:,:),1);

            end

        % Total power 
        elseif contains(metric,'tp')

            % Compute total power 
            eeg_features = squeeze(sum(power,1));

        % Root mean square frequency
        elseif contains(metric,'rmsf')

            % Compute root mean square frequency
            eeg_features = squeeze(sqrt(sum(repmat...
                ((f_vector'.^2),1,n_pnts).*power,1)));
        end


        %---------------------------------------------------------    
        % Remove outliers 
        %---------------------------------------------------------

        % Remove samples that are 10 times greater
        % than the signal's standard deviation 
        
        % Reshape feature matrix to have one
        % feature at each column 
        eeg_features = reshape(eeg_features, ...
            [n_pnts, numel(eeg_features(1,:,:))]);
        
        % Go through columns and replace 
        % outliers at each column 
        for c = 1 : size(eeg_features,2)
            
            col = eeg_features(:,c);
            col(col>10*std(col)) = 10*std(col);
            eeg_features(:,c) = col;
         
        end

        % Reshape feature matrix into its 
        % original dimension 
        eeg_features = reshape(eeg_features, ...
            [n_pnts,n_chans,n_bands]);

        %---------------------------------------------------------    
        % Mirror padd the signal before convolution  
        %---------------------------------------------------------

        % Assign first and last indices
        % according to the current sampling
        % frequency of the data 
        if strcmp(tf_method,'welch')
            f = first; l = last;
        elseif strcmp(tf_method,'wavelet')
            f = first_eeg; l = last_eeg;
        end

        % Mirror padd the signal at a length
        % equal to kernal size + 1 (seconds)
        eeg_features = eeg_features(f:l,:,:);

        % Padd feature in the pre-direction 
        padsize = max(f-1,kernel_size*fs_eeg);
        eeg_features = padarray(eeg_features, ...
            padsize,'symmetric','pre');

        % Padd feature in the post-direction
        padsize = max(n_pnts-l,kernel_size*fs_eeg);
        eeg_features = padarray(eeg_features, ...
            padsize,'symmetric','post');

        %---------------------------------------------------------    
        % Convolution with HRF or delay
        %---------------------------------------------------------

        switch eeg_shift 

            case 'conv'

                plotting_shift = 'convolved';
                eeg_features_delay = convolve_features(eeg_features, ...
                    fs_eeg,delays,kernel_size);

                % Cut the eeg predictor matrices to   
                % match the bold acquisition times 
                eeg_features = eeg_features(f:l,:,:);
                eeg_features_delay = eeg_features_delay(f:l,:,:,:);

            case 'delay'

                plotting_shift = 'delayed';
                eeg_features_delay = delay_features(eeg_features, ...
                    fs_eeg,delays,[f,l]); 

                % Cut the eeg predictor matrix to   
                % match the bold acquisition times 
                eeg_features = eeg_features(f:l,:,:);

            case ''

                % Cut the eeg predictor matrices to   
                % match the bold acquisition times 
                eeg_features = eeg_features(f:l,:,:);

        end  
        
        n_pnts_eeg = length(f:l);

        %---------------------------------------------------------    
        % Prune according to the task design (case sub_task)
        %---------------------------------------------------------
        if sub_task ~= ""
            
            sub_task_design = sub_task_design(f:l);
            eeg_features = eeg_features(logical(sub_task_design),:,:);
            
            if ~strcmp(eeg_shift,"")
                
                eeg_features_delay = eeg_features_delay...
                    (logical(sub_task_design),:,:,:);
                
            end
            
            n_pnts_eeg = size(eeg_features,1);
        end
        

        %---------------------------------------------------------    
        % Downsample to intermediate frequency
        %---------------------------------------------------------

        % Keep a version of the features in  
        % the original EEG sampling frequency 
        eeg_features_eeg_fs = eeg_features;
        
        switch tf_method

            case 'welch'

               % Assign time vector and # of points 
               n_pnts = length(first:last);
               time = 0 : 1/fs_eeg : (n_pnts-1)/fs_eeg;

            case 'wavelet'

                % Assign time vectors and # of points 
                time_eeg = 0 : 1/fs_eeg : (n_pnts_eeg-1)/fs_eeg;      
                time = 0 : 1/fs : (n_pnts_eeg-1)/fs_eeg; 
                n_pnts = length(time);

                % Use spline interpolation to downsample features 
                eeg_features = permute(eeg_features, [3 2 1]);         
                eeg_features = spline(time_eeg, eeg_features, time);         
                eeg_features = permute(eeg_features, [3 2 1]);                

                if ~strcmp(eeg_shift,"")
                    
                    eeg_features_delay = ...
                        permute(eeg_features_delay, [4 3 2 1]);    
                    eeg_features_delay = ...
                        spline(time_eeg, eeg_features_delay, time);
                    eeg_features_delay = ...
                        permute(eeg_features_delay, [4 3 2 1]);  
                    
                end
                
        end 


        %---------------------------------------------------------    
        % Normalize features (0 mean, 1 std)
        %---------------------------------------------------------

        % Normalize every feature so as to have
        % zero mean and standard deviation one
        eeg_features_norm = zscore(eeg_features);

        if ~strcmp(eeg_shift,"")
       
            eeg_features_delay_norm = ...
                zscore(reshape(eeg_features_delay, ...
                [n_pnts, numel(eeg_features_delay(1,:,:,:))]));
            eeg_features_delay_norm = ...
                reshape(eeg_features_delay_norm, ...
                size(eeg_features_delay));
            
        end

        %---------------------------------------------------------    
        % Plot features (after normalization)
        %---------------------------------------------------------
        
        if ~strcmp(eeg_shift,"") && flag.report ~= 0
            
            for c = 1 : n_chans

                for b = 1 : n_bands

                    figure_name = strcat(id_bands(b), ...
                        ' EEG features of'," ",subject, ...
                        ', channel'," ", num2str(c,'%01d')); 
                    figure('Name',figure_name)

                    plot(time, eeg_features_norm(:,c,b)); hold on; 
                    plot(time, eeg_features_delay_norm(:,c, ...
                        plotting_delay,b));

                    title(strcat(id_bands(b), ' feature of'," ", ...
                        subject,'channel'," ", num2str(c,'%01d'), ...
                        ' (demeaned, normalized)')); 
                    xlabel('Time(s)'); ylabel('Power (\muV^2)');
                    legend(strcat(id_bands(b),' feature'), ...
                        strcat(id_bands(b),' feature,', ...
                        plotting_shift,', delay'," " , ...
                        num2str(delays(plotting_delay)), 's'));

                    img_out = strcat(id_bands(b),num2str(delays ...
                        (plotting_delay)),'sec',num2str(c),'chan.png');
                    saveas(gcf,fullfile(path_img_out(s),img_out));

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
            [n_pnts_eeg,numel(eeg_features_eeg_fs(1,:,:))]);
        
        if ~strcmp(eeg_shift,"")
            
            eeg_features_delay_norm = ...
                reshape(eeg_features_delay_norm, ...
                [n_pnts,numel(eeg_features_delay_norm(1,:,:,:))]);
            
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
        dlmwrite(fullfile(path_data_out(s), ...
            feature_out),eeg_features_norm);

        dlmwrite(fullfile(path_data_out(s), ...
            feature_out_eeg_fs),eeg_features_eeg_fs);
        
        if ~strcmp(eeg_shift,"")
            dlmwrite(fullfile(path_data_out(s), ...
                feature_out_delay),eeg_features_delay_norm);
        end

    end % finish looping through subjects 

end % finish looping through metrics 