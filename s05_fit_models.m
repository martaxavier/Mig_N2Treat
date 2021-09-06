% This script calls either kfold_cv, kfold_cv_nondep, kfold_cv_blocked,
% k_fold_cv_sessions 
        
%--------------------------------------------------------    
% Go through metrics
%-------------------------------------------------------- 

for m = 1 : length(metrics)  

   metric = metrics(m);

   % Get parameters for 
   % current metric 
   get_metric_pars;

    % If cross-validation is to be performed
    % across sessions, it won't be necessary 
    % to compute the optimal cv parameters 
    if ~strcmp(cv_method, 'sessions') && ...
        ~strcmp(cv_method, 'one_class')
    
        % Load estimated order of the ACF for current metric 
        if strcmp(bold_shift, 'deconv')
            acf_order_in = strcat('acf_order_deconv_', sessions, '.mat');
        else
            acf_order_in = strcat('acf_order_', sessions', '.mat');
        end
        load(fullfile(path_pars, acf_order_in(se)), 'acf_order');
        
    end 
               
    % Define input EEG and BOLD data, according to
    % the current metric
    eeg_in = strcat(eeg_metric, '_', 'eeg_feature', ...
        '_', eeg_shift, '.txt');
    bold_in = strcat('bold_preproc', '_', bold_shift, '.txt');

    if contains(eeg_in, '_.')
        eeg_in = replace(eeg_in, '_.', '.');
    end
    if contains(bold_in, '_.')
        bold_in = replace(bold_in, '_.', '.');
    end
            

   %---------------------------------------------------    
   % Go through subjects
   %----------------------------------------------------

   n_subjects = length(subjects);
   n_sessions = length(sessions);
   
   if strcmp(cv_method, 'one_class')
             
        % Create outpur directory, if not existent
        if ~exist(path_data_out, 'dir')
           mkdir(path_data_out);
        end  
        
        if ~exist(path_hc, 'dir')
            mkdir(path_hc);
        end
            
        eeg = dlmread(fullfile(path_eeg_in(1, 1), eeg_in));
        bold = dlmread(fullfile(path_bold_in(1, 1), bold_in));

        eeg = zeros([size(eeg) n_subjects*n_sessions]);
        bold = zeros([length(bold) n_subjects*n_sessions - 1]); 
        s_se = 1;
        
   end
   
   for s = 1 : n_subjects

        subject = subjects(s);          

        % If cross-validation is to be performed
        % across sessions, it won't be necessary 
        % to compute the optimal cv parameters 
        if strcmp(cv_method, 'sessions')
            
            % Create outpur directory, if not existent
            if ~exist(path_data_out(s), 'dir')
               mkdir(path_data_out(s));
            end 
            
            eeg = dlmread(fullfile(path_eeg_in(s, 1), eeg_in));
            bold = dlmread(fullfile(path_bold_in(s, 1), bold_in));
            
            eeg = cat(3, eeg, zeros([size(eeg) 2]));
            bold = cat(2, bold, zeros([length(bold) 2]));
            
            % Go through sessions
            for se = 2 : n_sessions
                
                % Load input EEG and BOLD data 
                eeg(:, :, se) = dlmread(fullfile(path_eeg_in(s, se), eeg_in));
                bold(:, se) = dlmread(fullfile(path_bold_in(s, se), bold_in));
                
            end % sessions 
            
        elseif strcmp(cv_method, 'one_class')
            
            % Go through sessions 
            for se = 1 : n_sessions
                
                EEG = dlmread(fullfile(path_eeg_in(s, se), eeg_in));
                eeg(:, :, s_se) =  dlmread(fullfile(path_eeg_in(s, se), eeg_in));
                bold(:, s_se) = dlmread(fullfile(path_bold_in(s, se), bold_in));
                s_se = s_se + 1; 
                
            end % sessions 
            
        else
            
            % Create output directory, if not existent 
            if ~exist(path_data_out(s, se), 'dir')
               mkdir(path_data_out(s, se));
            end 

            % Load input EEG and BOLD data 
            eeg = dlmread(fullfile(path_eeg_in(s, se), eeg_in));
            bold = dlmread(fullfile(path_bold_in(s, se), bold_in));
            
        end

        % Perform cross-validation to fit the current model
        % Use the cv method explicit in the variable cv 
        switch cv_method

            case 'regular'
   
                display(strcat('Fitting model for', " ", ...
                    subject, ',', " ", metric, ',', " ", ...
                    cv_method, ',', " ", reg_models, '...'));
        
                [model, optimal] = kfold_cv_par_v2(eeg, bold,...
                    'k', 10, 'val2learn', 0.2, 'regress', reg_models);

            case 'nondep'
                
                 display(strcat('Fitting model for', " ", ...
                    subject, ',', " ", metric, ',', " ", ...
                    cv_method, ',', " ", reg_models, '...'));               
                
                idx = find(ismember(table2array(acf_order), subject));
                order = table2array(acf_order(idx, 2));   
                [model, optimal] = kfold_cv_nondep_par_v3...
                    (eeg, bold, 'k', 10, 'v', 0.3, 'regress', reg_models, ...
                    'autocorr', order, 'sizx', [length(bold) prod(dim)]);

            case 'blocked'
                
                display(strcat('Fitting model for', " ", ...
                    subject, ',', " ", metric, ',', " ", ...
                    cv_method, ',', " ", reg_models, '...'));                

                [model, optimal] =  kfold_cv_blocked_par_v3(eeg, bold, ...
                    path_pars, 'k', 10, 'v', 0.3, 'regress', reg_models, ...
                    'sizx', [length(bold) prod(dim)]);
                
            case 'sessions'
                
                display(strcat('Fitting model for', " ", ...
                    subject, ',', " ", metric, ',', " ", ...
                    cv_method, ',', " ", reg_models, '...'));                
                   
                [optimal] = kfold_cv_sessions_par(eeg, bold, path_pars, ...
                    'regress', reg_models, 'sizx', [size(bold, 1) ...
                    prod(dim)]);

        end 

        % Define output data according to current metric 
        model_out = strcat(metric, '_', 'model_', cv_method, '.mat'); 
        model_folds_out = strcat(metric, '_',  'model_folds_', cv_method, '.mat'); 

        % Save model
        if ~strcmp(cv_method, 'sessions') && ~strcmp(cv_method, 'one_class')
            save(fullfile(path_data_out(s), model_out), 'model');
            save(fullfile(path_data_out(s), model_folds_out), 'optimal');            
        elseif ~strcmp(cv_method,  'one_class')
            save(fullfile(path_data_out(s), model_folds_out), 'optimal');
        end

   end % finish looping through subjects 
   
   switch cv_method
       
       case 'one_class'
           
        display(strcat('Fitting model for', " ", metric, ',', ...
            " ", cv_method, ',', " ", reg_models, '...'));           
           
         [ses_coef, clusters, optimal] = kfold_one_class_par(eeg, bold, ...
             path_pars, path_hc, 'hc_distance', hc_distance, 'v', 0.3, 'n', 4); 
         
         % Save model
         model_one_class_out = strcat(metric, '_', 'model_', cv_method, '.mat');  
         ses_coef_out = strcat(metric, '_', 'ses_coef_', cv_method, '.mat'); 
         clusters_out = strcat(metric, '_', 'clusters_', cv_method, '.mat'); 
         save(fullfile(path_data_out, model_one_class_out), 'optimal');
         save(fullfile(path_data_out, ses_coef_out), 'ses_coef');
         save(fullfile(path_data_out, clusters_out), 'clusters');
           
   end

end % finish looping through metrics
