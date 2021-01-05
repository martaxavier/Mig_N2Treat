% This script calls either kfold_cv, kfold_cv_nondep, kfold_cv_blocked      
        
%--------------------------------------------------------    
% Go through metrics
%-------------------------------------------------------- 

for m = 1 : length(metrics)  

   metric = metrics(m);

   % Get parameters for 
   % current metric 
   get_metric_pars;
   
%     % Continue if metric is not yet supported for model fitting 
%     if strcmp(reg_model,'l21_1') && n_bands == 1
%         disp(char(strcat('Model wont be fitted for metric', ...
%             {' '}, metric,' because this metric is not', ...
%             ' supported for regression method', {' '}, ...
%             reg_model, '...')));
%         continue
%     end
    
    % Load optimal CV parameters K and V for current metric 
    cv_pars_in = strcat(reg_model,'_',cv_method,'_',metric,'.mat');
    load(char(fullfile(path_pars_in,cv_pars_in)),'cv_pars');
    k = cv_pars.K; v = cv_pars.V;

    % Load estimated order of the ACF for current metric 
    if strcmp(bold_shift,'deconv')
        acf_order_in = 'acf_order_deconv.mat';
    else
        acf_order_in = 'acf_order.mat';
    end
    load(char(fullfile(path_pars_in, ...
        acf_order_in)),'acf_order');

   %---------------------------------------------------    
   % Go through subjects
   %----------------------------------------------------

   for s = 1 : length(subjects)

        subject = subjects(s);          

        disp(char(strcat('Fitting model for', {' '}, ...
            subject, ',', {' '}, metric, ',', {' '}, ...
            cv_method, ',', {' '}, reg_model, '...')));

        % Create output directory, if not existent 
        if ~exist(char(path_data_out(s)), 'dir')
           mkdir(char(path_data_out(s)));
        end 

        % Define input EEG and BOLD data, according
        % to current metric
        eeg_in = strcat(eeg_metric,'_', ...
            'eeg_feature','_',eeg_shift,'.txt');
        bold_in = strcat('bold_processed', ...
            '_',bold_shift,'.txt');

        if contains(char(eeg_in),'_.')
            eeg_in = replace(eeg_in,'_.','.');
        end
        if contains(char(bold_in),'_.')
            bold_in = replace(bold_in,'_.','.');
        end

        % Load input EEG and BOLD data 
        eeg=dlmread(char(fullfile(path_eeg_in(s),eeg_in)));
        bold=dlmread(char(fullfile(path_bold_in(s),bold_in)));

        % Perform crossvvalidation to fit the current model
        % Use the cv method explicit in the variable cv 
        switch cv_method

            case 'regular'

                [model,optimal] = kfold_cv_par_v2(eeg,bold,...
                    'k',k,'v',v,'regress',reg_model);

            case 'nondep'
                
                idx = find(ismember(table2array(acf_order),subject));
                order = table2array(acf_order(idx,2));
                [model,optimal] = kfold_cv_nondep_par_v3...
                    (eeg,bold,'k',k,'regress',reg_model, ...
                    'autocorr',order,'sizx', ...
                    [length(bold) prod(dim(1:2)) dim(3)]);

            case 'blocked'

                idx = find(ismember(table2array(acf_order),subject));
                order = table2array(acf_order(idx,2));
                [model,optimal] =  kfold_cv_blocked_par_v3(eeg,bold,...
                    'k',k,'val2learn',v,'regress',reg_model,...
                    'autocorr',order,'sizx', ...
                    [length(bold) prod(dim(1:2)) dim(3)]);

        end 

        % Define output data according to current metric 
        model_out=strcat(metric,'_','model_', ...
           cv_method,'.mat'); 
        model_folds_out=strcat(metric,'_', ...
            'model_folds_',cv_method,'.mat'); 

        % Save model and optimal structs in .mat files 
        save(char(fullfile(path_data_out(s),model_out)),'model');
        save(char(fullfile(path_data_out(s),model_folds_out)),'optimal');

   end % finish looping through subjects 

end % finish looping through metrics
