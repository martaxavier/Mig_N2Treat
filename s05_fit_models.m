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
    if ~strcmp(cv_method, 'sessions')
        
        % Load optimal CV parameters K and 
        % V for current metric 
        cv_pars_in = strcat(reg_models, '_', ...
            cv_method, '_', metric, '.mat');
        load(fullfile(path_pars, cv_pars_in), ...
            'cv_pars');
        k = cv_pars.K; v = cv_pars.V;

        % Load estimated order of the ACF for current metric 
        if strcmp(bold_shift, 'deconv')
            acf_order_in = strcat('acf_order_deconv_', sessions, '.mat');
        else
            acf_order_in = strcat('acf_order_', sessions', '.mat');
        end
        load(fullfile(path_pars, ...
            acf_order_in), 'acf_order');
        
    end 

   %---------------------------------------------------    
   % Go through subjects
   %----------------------------------------------------

   for s = 1 : length(subjects)

        subject = subjects(s);          

        display(strcat('Fitting model for', " ", ...
            subject, ',', " ", metric, ',', " ", ...
            cv_method, ',', " ", reg_models, '...'));
              
        % Define input EEG and BOLD data, according
        % to current metric
        eeg_in = strcat(eeg_metric, '_', ...
            'eeg_feature', '_', eeg_shift, '.txt');
        bold_in = strcat('bold_preproc', ...
            '_', bold_shift, '.txt');

        if contains(eeg_in, '_.')
            eeg_in = replace(eeg_in, '_.', '.');
        end
        if contains(bold_in, '_.')
            bold_in = replace(bold_in, '_.', '.');
        end
        
        % If cross-validation is to be performed
        % across sessions, it won't be necessary 
        % to compute the optimal cv parameters 
        if ~strcmp(cv_method, 'sessions')

            % Create output directory, if not existent 
            if ~exist(path_data_out(s, se), 'dir')
               mkdir(path_data_out(s, se));
            end 

            % Load input EEG and BOLD data 
            eeg = dlmread(fullfile(path_eeg_in(s, se), eeg_in));
            bold = dlmread(fullfile(path_bold_in(s, se), bold_in));
            
        else
            
            % Create outpur directoru, if not existent
            if ~exist(path_data_out(s), 'dir')
               mkdir(path_data_out(s));
            end 
            
            eeg = dlmread(fullfile(path_eeg_in(s, 1), eeg_in));
            bold = dlmread(fullfile(path_bold_in(s, 1), bold_in));
            
            eeg = cat(3, eeg, zeros([size(eeg) 2]));
            bold = cat(2, bold, zeros([length(bold) 1]));
            
            % Go through sessions
            for se = 2 : length(sessions)
                
                % Load input EEG and BOLD data 
                eeg(:, :, se) = dlmread(fullfile(path_eeg_in(s, se), eeg_in));
                bold(:, se) = dlmread(fullfile(path_bold_in(s, se), bold_in));
                
            end % sessions 
            
        end

        % Perform cross-validation to fit the current model
        % Use the cv method explicit in the variable cv 
        switch cv_method

            case 'regular'

                [model, optimal] = kfold_cv_par_v2(eeg, bold,...
                    'k', k, 'val2learn', v, 'regress', reg_models);

            case 'nondep'
                
                idx = find(ismember(table2array(acf_order), subject));
                order = table2array(acf_order(idx, 2));   
                [model, optimal] = kfold_cv_nondep_par_v3...
                    (eeg, bold, 'k', k, 'regress', reg_models, ...
                    'autocorr', order, 'sizx', ...
                    [length(bold) prod(dim)]);

            case 'blocked'

                idx = find(ismember(table2array(acf_order), subject));
                order = table2array(acf_order(idx, 2));
                [model, optimal] =  kfold_cv_blocked_par_v3(eeg, bold, ...
                    'k', 6, 'v', 6, 'regress', reg_models,...
                    'autocorr', order, 'sizx', ...
                    [length(bold) prod(dim)]);
                
            case 'sessions'
                   
                [optimal] = kfold_cv_sessions_par(eeg, bold, ...
                    'regress', reg_models, 'sizx', [size(bold, 1) ...
                    prod(dim)]);

        end 

        % Define output data according to current metric 
        model_out = strcat(metric, '_', 'model_', ...
           cv_method, '.mat'); 
        model_folds_out = strcat(metric, '_', ...
            'model_folds_', cv_method, '.mat'); 

        % Save model and optimal structs in .mat files 
        if ~strcmp(cv_method, 'sessions')
            save(fullfile(path_data_out(s), model_out), 'model');
        end
        save(fullfile(path_data_out(s), model_folds_out), 'optimal');

   end % finish looping through subjects 

end % finish looping through metrics
