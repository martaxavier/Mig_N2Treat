% ------------------------------------------------------------
% Go through metrics
% ------------------------------------------------------------
        
for m = 1 : length(metrics)
    
    metric = metrics(m);
    get_metric_pars;
   
    if strcmp(reg_models,'l21_1') && n_bands == 1
        disp(strcat('CV pars wont be estimated for metric', ...
            " ", metric, ' because this metric is not', ...
            ' supported for regression model', " ", ...
            reg_models, '...'));
        continue
    end
    
    disp(strcat('Estimating CV pars for regression model', ...
        " ", reg_models, ', metric', " ", metric, '...'));
    
    % Define input ACF data 
    if strcmp(bold_shift, 'deconv')
        acf_order_in = 'acf_order_deconv.mat';
    else
        acf_order_in = 'acf_order.mat';
    end
    
    % Load order of the ACF model for the current subject 
    load(fullfile(path_pars_in, acf_order_in), 'acf_order');

    % Pre-allocate matrices 
    order = zeros(length(subjects), 1);
    X = cell(length(subjects), 1);
    Y = cell(length(subjects), 1);
    
    % ------------------------------------------------------------
    % Go through subjects and get X and Y
    % ------------------------------------------------------------
    for s = 1 : length(subjects)

        subject = subjects(s);
        
        % Get order of the ACF for current subject and store it
        idx = find(ismember(table2array(acf_order), subject));
        order(s) = table2array(acf_order(idx, 2));
    
        % Define input EEG and BOLD data, according
        % to current metric
        eeg_in = strcat(eeg_metric, '_', ...
            'eeg_feature', '_', eeg_shift, '.txt');
        bold_in = strcat('bold_processed', ...
            '_', bold_shift, '.txt');

        if contains(eeg_in,'_.')
            eeg_in = replace(eeg_in, '_.', '.');
        end
        if contains(bold_in,'_.')
            bold_in = replace(bold_in, '_.', '.');
        end

        % Get BOLD signal and EEG features for current subject
        X{s} = dlmread(char(fullfile(path_eeg_in(s), eeg_in)));
        Y{s} = dlmread(char(fullfile(path_bold_in(s), bold_in)));

    end
  
    siz_X = [size(X{s},1) prod(dim(1:2)) dim(3)];

    % ------------------------------------------------------------
    % Get stats and best k-v pair for all subjects 
    % ------------------------------------------------------------

    [best_pair, stats] = get_cv_pars(X, Y, cv_method,...
        reg_models, order, siz_X);

    % ------------------------------------------------------------
    % Save output matrices with optimal parameters  
    % ------------------------------------------------------------

    % Save best K-V pair in outputstruct
    cv_pars.K = best_pair(1); cv_pars.V = best_pair(2);

    % Save matrix of best parameters for this cv method  
    cv_pars_out = strcat(reg_models, '_', cv_method, '_', metric, '.mat');
    save(fullfile(path_pars, cv_pars_out), 'cv_pars');    

end

%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% [best_pair,stats] = get_cv_pars(X, Y, cv, reg_models, order, siz_X)       
% ============================================================

function [best_pair,stats] = get_cv_pars(X, Y, cv, reg_models, order, siz_X)

%   Input parameters:
%
%     X                 EEG feature of each subject
%     Y                 BOLD of each subject
%     cv                The CV method    
%     regress           The regression method 
%     autocorr          The order of the ACF
%     sizx              The size of the problem
%
%   Outputs:
%
%       bestpairs       The optimal K-V pair. K is the total number of
%                       folds, V is the validation to learn set ratio
%
%
%       stats            The struct that contains information about the
%                        sequence of model fits corresponding to the
%                        the searched parameters. stats contains the
%                        following fields:
%  
%           nmse         Normalized mean squared error of the 
%                        fits for each of the combinations tested 
%
%           df           Degrees of freedom of the fits for 
%                        each of the combinations tested 
%
%           bic          Bayesian information criterion of the 
%                        fits for each of the combinations tested 

tic 

n_subjects = length(X);

% ------------------------------------------------------------
% Choose the K and V values to be tested 
% ------------------------------------------------------------

% Typical values of K or V
% range from 5 to 15, however
% because of the non-dependent 
% CV, a partition of 5 would 
% result in too big test set 
K = 15:5:15; 
n_k = length(K);

V = 10:5:15; 
n_v = length(V);

% Linear indices 
idxs = 1:n_k*n_v;
n_pairs = length(idxs);

% Subscript indices 
[k_idxs, v_idxs] = ...
    ind2sub([n_k,n_v], idxs);

nmse = zeros(n_subjects, n_pairs); 
df = nmse; bic = nmse;

% ------------------------------------------------------------
% Go through K-V pairs 
% ------------------------------------------------------------
for p = 1 : n_pairs
    
    k = K(k_idxs(p));
    v = V(v_idxs(p));
        
    % ------------------------------------------------------------
    % Go through  subjects 
    % ------------------------------------------------------------
    for s = 1 : n_subjects
        
        disp(strcat('Subject', " ", num2str(s), '...'));
      
        EEG = X{s};
        BOLD = Y{s};
        
        % Center and scale EEG
        [EEG,~,~] = zscore(EEG, 1);

        % Center BOLD
        muBOLD = mean(BOLD);
        BOLD = bsxfun(@minus, BOLD, muBOLD);
        
        switch cv 
            
            case 'regular'
                
                [model,~] = kfold_cv_par_v2(EEG, BOLD, ...
                    'k', k, 'v', v, 'regress', reg_models); 
                
            case 'nondep'
                
                [model,~] = kfold_cv_nondep_par_v3...
                    (EEG, BOLD, 'k', k, 'v', v, ...
                    'regress', reg_models, ...
                    'autocorr',order(s),'sizx',siz_X);   
                
            case 'blocked'
                
                [model,~] = kfold_cv_blocked_par_v2...
                    (EEG, BOLD, 'k', k, 'v', v, ...
                    'regress', reg_models, 'autocorr', ...
                    order(s), 'sizx', siz_X);    
                
        end 
    
        nmse(s,p) = model.nmse;
        df(s,p) = model.df;
        bic(s,p) = model.bic;  
        
    end % finish looping through subjects 

end % finish looping through k-v pairs 

% ------------------------------------------------------------
% Interpolate remaining K-V pairs
% ------------------------------------------------------------

% Average model results across subject dimension
nmse_avg = mean(nmse, 1);
df_avg = mean(df, 1); 
bic_avg = mean(bic, 1);

% Reshape results into matrices for plotting 
nmse_mat = reshape(nmse_avg, [n_k, n_v]);
df_mat = reshape(df_avg, [n_k, n_v]);
bic_mat = reshape(bic_avg, [n_k, n_v]);

% Plot surf of nmse, df and bic 
[K_2d, V_2d] = meshgrid(K, V);

figure;surf(K_2d, V_2d, nmse_mat);
title('NMSE'); xlabel('K'); ylabel('V');

figure;surf(K_2d, V_2d, df_mat);
title('DF'); xlabel('K'); ylabel('V');

figure;surf(K_2d, V_2d, bic_mat);
title('BIC'); xlabel('K'); ylabel('V');

% ------------------------------------------------------------
% Find best K-V pair 
% ------------------------------------------------------------

% Sort indices according to bic 
[~,sortedindices] = sort(bic_avg);
best_pair = [K(k_idxs((sortedindices(1)))), ...
    V(v_idxs(sortedindices(1)))];

% ------------------------------------------------------------
% Prepare outputs 
% ------------------------------------------------------------

% Prepare stats struct 
stats.k = K(k_idxs((sortedindices)))';
stats.v = V(v_idxs((sortedindices)))';
stats.nmse = nmse_avg(sortedindices)';
stats.df = df_avg(sortedindices)';
stats.bic = bic_avg(sortedindices)';
stats.summary = struct2table(stats);

toc

end 