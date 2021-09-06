function [optimal] = kfold_cv_sessions_par(EEG, BOLD, path_pars, varargin)

%   [model, optimal] = kfold_cv_sessopns_par(EEG, BOLD, ...) performs,  
%   nested k-fold cross-validation, across sessions, to obtain 
%   the best model for the input EEG and BOLD data
%
%   Input data:
%
%     EEG              The model's X 
%     BOLD             The model's Y
%   
%   Optional input parameters:  
%
%     'v'              The validation to learn fraction of the 
%                      inner CV loop 
%     'n'              The number of cyles in the inner CV loop
%     'regress'        The regression method used to fit the model
%     'rho'            The complexity parameter rho of the L21+1/L2+1 fit 
%     'lambda'         The complexity parameter lambda of the L21+1/L2+1 fit
%     'sizx'           The size of the problem
%   
%   Outputs:
%
%     model            model is a struct with information corresponding to 
%                      the best fit found, and contains the following
%                      fields:
%
%       'rho'          The regularization parameter rho of the L21+1 fit
%       'lambda'       The complexity parameter lambda of the L21+1 fit 
%       'efp'          The activation pattern of the final fit
%       'yhat'         The Y hat corresponding to the final fit 
%       'df'           The number of degrees of freedom of the final fit 
%       'nmse'         The normalized mean squared error of the final fit 
%       'bic'          The bayesian inference criterion of the final fit 
%

tic 

% ------------------------------------------------------------ 
% Sanity check and process input parameters 
% ------------------------------------------------------------ 

% X a real 3D matrix
if length(size(EEG)) ~= 3 || ~isreal(EEG)
    error('EEG is not a real 3D matrix');
end

% If Y is a #sessions-row vector, 
% convert to a #sessions-column vector
if size(BOLD, 1) < size(BOLD, 2)
    BOLD = BOLD';
end

% Y a vector, same length as the columns of X
if ~ismatrix(BOLD) || ~isreal(BOLD) || length(size(BOLD)) ~= 2
    error('BOLD is not a real 2D matrix');
end

% Number of regressors in the model
n_features = size(EEG, 2);

% Number of sessions = outer folds
K = size(EEG, 3);

% ------------------------------------------------------------ 
% Sanity check and process optional parameters 
% ------------------------------------------------------------ 

% Assign default values for each optional parameter
pnames = {'v' 'n' 'regress' 'rho' 'lambda' 'sizx'}; 
dflts  = {0.2 4 'l2_1' [] [] []};

% Assign variables corresponding to optional parameters 
[V, N, method, rho, lambda, siz_X] = ...
    internal.stats.parseArgs(pnames, dflts, varargin{:});

% Check if user supplied lambda and rho
user_supplied_rho = true;
user_supplied_lambda = true;

if isempty(rho)
    user_supplied_rho = false;
end    

if isempty(lambda)
    user_supplied_lambda = false;
end

if ~isequal(user_supplied_rho,...
        user_supplied_lambda)
    error('Must supply both rho and lambda, or neither');
end
        
% ------------------------------------------------------------ 
% Nested cv procedure  
% ------------------------------------------------------------ 

% Allocate arrays of model
% parameters and performance
% for each test/train pair
opt_lambda =	zeros(K, 1); 
opt_rho =       opt_lambda;
opt_coef =      zeros(n_features+1, K);
opt_df =        opt_lambda; 

opt_bic_train =  opt_lambda; 
opt_mse_train =  opt_lambda;
opt_nmse_train = opt_lambda; 
opt_corr_train = opt_lambda; 

opt_bic_test =  opt_lambda; 
opt_mse_test =  opt_lambda;
opt_nmse_test = opt_lambda; 
opt_corr_test = opt_lambda; 

% In case the user hasn't supplied rho and lambda (must 
% supply both or neither), fix rho and retreive pre-established
% lambda range of interest
if ~user_supplied_rho   
    
    if strcmp(method, 'l2_1') || strcmp(method, 'l21_1')
    
        % Rho or alpha parameter
        rho = 0.6;

        % Retreive pre-established lambda range of interest
        load(fullfile(path_pars, strcat(method, ...
            '_Lambdas.mat')), 'lambda');
        
    elseif strcmp(method, 'rf')
        
        % Number of trees
        rho = 100;
        
        % Retreive pre-established range of number of features 
        % to be used in each split
        n_features_split_def = round(n_features/3);
        llim = n_features_split_def - round(n_features/20);
        hlim = n_features_split_def + round(n_features/20);
        lambda = round(linspace(llim, hlim, 5));
        
    end
    
end

% Number of lambda values 
n_lambda = length(lambda);

%-------------- Begin outer loop ---------------%

for k = 1 : K

    % Assign train set indices
    ind_train = k;
    
    % Assign test set indices
    ind_test = find(1 : K ~= k);

    % Assign train and test X (EEG) and Y (BOLD) variables 
    X_train = EEG(:, :, ind_train); y_train = BOLD(:, ind_train);
    X_test = EEG(:, :, ind_test); y_test = BOLD(:, ind_test);
    
    % Linearize test sets 
    X_test = reshape(permute(X_test, [1 3 2]), ...
        [size(X_test, 1)*size(X_test, 3) size(X_test, 2)]);
    y_test = y_test(:);
    
    % Obtain train and test set sizes
    siz_train = size(y_train, 1);
    siz_test = size(y_test, 1);
    
    % Allocate bic and mse matrices for the learn and 
    % for the val set, each inner iteration through cols
    bic_val = zeros(n_lambda, N); nmse_val = bic_val;
    df_inner = zeros(n_lambda, N);
    oob_error = zeros(n_lambda, N);

    %-------------- Begin inner loop --------------%
        
    % The inner loop 
    % has N iterations 
    parfor n = 1 : N
        
        % Assign broadcast variables to loop variables for efficiency 
        % Large broadcast variables can cause significant communication 
        % between client and workers and increase parallel overhead
        % replace them for temporary variables, created inside the loop 
        EEG_par = EEG; BOLD_par = BOLD; lambda_par = lambda;        
        
        % Circular holdout assignment 
        ind_start = randi([1 siz_train]);
        indices_in = sort(crossvalind('Holdout', siz_train, V));
    
        % Assign learning set indices 
        ind_learn = zeros(size(indices_in));
        ind_learn(ind_start : end) = indices_in(1 : size(ind_learn) - ind_start + 1);
        ind_learn(1 : ind_start - 1) = indices_in(size(ind_learn) - ind_start + 2 : size(ind_learn));        
        
        % Assign validation set indices 
        ind_val = (~ind_learn);
        ind_val = find(ind_val);
        siz_val = length(ind_val);
        
        ind_learn = find(ind_learn);         
        
        % Assign learning and validation variables 
        X_learn = squeeze(EEG_par(ind_learn, :, ind_train));
        y_learn = squeeze(BOLD_par(ind_learn, ind_train));
        X_val = squeeze(EEG_par(ind_val, :, ind_train)); 
        y_val = squeeze(BOLD_par(ind_val, ind_train));
            
        % Screen the input method
        if strcmp(method, 'l21_1')

           % L21+1 fit
            [betas_par, stats_par] = regress_L21_1(X_learn, y_learn, ...
                siz_X, 'Rho', rho, 'Lambda', lambda_par);
            
            [~, col] = find(betas_par); 
            df = accumarray(col, 1); 
            df(setdiff(1:n_lambda, col))= 0; 
            df = flip(df);
            betas_par = flip(betas_par, 2);          

        elseif strcmp(method, 'l2_1')

            % L2+1 elastic-net regression fit 
            [betas_par, stats_par] = lasso(X_learn, y_learn, ...
                'Alpha', rho, 'Lambda', lambda_par);
            
            df = stats_par.DF; df = flip(df)'; 
            betas_par = flip(betas_par, 2);                
            
        elseif strcmp(method, 'rf')
            
            oob_error_par = zeros(n_lambda, 1);

            for i = 1 : n_lambda
                
                % Rho is the # of trees and lambda the # of variables 
                % to select at random for each decision split
                rf_model_par = TreeBagger(rho, X_learn, y_learn, 'Method', ...
                    'regression', 'OOBPredictorImportance', 'on', ...
                    'NumPredictorsToSample', lambda_par(i));     
                oob_error_tmp = oobError(rf_model_par);
                oob_error_par(i) = oob_error_tmp(end);
            
            end

            oob_error(:, n) = oob_error_par;

        end
        
          
        if strcmp(method, 'l2_1') || strcmp(method, 'l21_1')
            
            % Save intercept values for current rho-lambda
            intercept = stats_par.Intercept; intercept=flip(intercept);
            y_hat_val = intercept + X_val*betas_par; 

            % Compute bic values for all rho-lambda 
            % pairs in the val and learn set     
            mse_val = sum((y_hat_val - y_val).^2)';     
            nmse_val(:, n) = mse_val/sum((y_val - mean(y_val)).^2);     
            bic_val(:, n) = log(siz_val).*df + ...
                siz_val.*log(mse_val ./ siz_val);

            df_inner(:, n) = df;
            
        end
        
    end
    
    %--------------- End inner loop ---------------%
    
    if strcmp(method, 'l2_1') || strcmp(method, 'l21_1')
    
        % Flag as ineligible, in 'trash_lam' (1) lambda values for which 
        % nmse values are, at any of the inner iterations, above a given 
        % inacceptable threshold, 'thresh_nmse' and (2) lambda values for 
        % which dof values are, at any of the inner iterations, zero
        % Flag, in 'trash_n', inner iterations in which the average nmse 
        % value in the validation set was above 'thresh_nmse' 
        thresh_nmse = 0.98; aux = nmse_val > thresh_nmse; 
        aux(:, mean(nmse_val, 1) > thresh_nmse - 0.05) = 0; 
        trash_lam = find(sum(aux, 2) >= 1); 
        [rows, ~] = find(df_inner == 0); 
        trash_lam = unique([trash_lam; unique(rows)]);

        % Average bic values without considering
        % iterations flagged as "trash" 
        aux = bic_val; aux(:, mean(nmse_val, 1) > thresh_nmse) = 0; 
        [rows,~,val] = find(aux);
        aux_val = accumarray(rows, val, [], @mean);

        % Find optimal rho-lambda for the current test/train pair,
        % i.e, find rho-lambda that minimizes sum of bics through  
        % all V inner iterations in the validation set
        aux_val(trash_lam) = Inf; [~, ind_opt] = min(aux_val); 
        if isempty(ind_opt); [~, ind_opt] = min(sum(bic_val, 2)); end
    
    elseif strcmp(method, 'rf')
        
        % Find optimal number of leafs that
        % resulted in minimum out-of-bag error
        [~, ind_opt] = min(sum(oob_error, 2)); 
        
    end
       
    % Save rho-lambda pair that minimizes
    % bic in learn and val sets combined 
    opt_lambda(k, 1) = lambda(ind_opt);
    opt_rho(k, 1) = rho;

    if strcmp(method, 'l21_1')

        % L21+1 fit for the current iteration of the outer CV loop 
        [betas, stats] = regress_L21_1(X_train, y_train, siz_X,...
        'Rho', opt_rho(k), 'Lambda', opt_lambda(k));
        betas(abs(betas) < 5e-4) = 0; 
        opt_df(k) = length(find(betas));
        
        % Compute the y hat for the test and training 
        % sets of the current outer iteration k
        y_hat_test = stats.Intercept + X_test*betas;
        y_hat_train = stats.Intercept + X_train*betas;  
        
        % Model coefficientrs of the current outer 
        % iteration k 
        opt_coef(:, k) = [stats.Intercept; betas];        

    elseif strcmp(method, 'l2_1')
        
        % L2+1 fit for the current iteration of the outer CV loop 
        [betas, stats] = lasso(X_train, y_train, 'Alpha', opt_rho(k), ...
            'Lambda', opt_lambda(k)); opt_df(k) = stats.DF; 
        
        % Compute the y hat for the test and training 
        % sets of the current outer iteration k
        y_hat_test = stats.Intercept + X_test*betas;
        y_hat_train = stats.Intercept + X_train*betas;
        
        % Model coefficientrs of the current outer 
        % iteration k 
        opt_coef(:, k) = [stats.Intercept; betas];   
        
    elseif strcmp(method, 'rf')
        
        rf_model = TreeBagger(opt_rho(k), X_train, y_train, ...
            'Method', 'regression', 'OOBPredictorImportance', 'on', ...
            'NumPredictorsToSample', opt_lambda(k)); 
        
        [y_hat_test, ~] = predict(rf_model, X_test); 
        [y_hat_train, ~] = predict(rf_model, X_train);  
    
    end
    
    % Model performance of the current outer 
    % iteration k, in the test and train sets 
    opt_mse_test(k) = sum((y_hat_test - y_test).^2);
    opt_mse_train(k) = sum((y_hat_train - y_train).^2);
    
    opt_bic_test(k) = log(siz_test).* opt_df(k) + ...
        siz_test.* log(opt_mse_test(k)./ siz_test);
    opt_bic_train(k) = log(siz_train).* opt_df(k) + ...
        siz_train.* log(opt_mse_train(k)./ siz_train);
    
    opt_nmse_test(k) = opt_mse_test(k)/ ...
        sum((y_test - mean(y_test)).^2);
    opt_nmse_train(k) = opt_mse_train(k)/ ...
        sum((y_train - mean(y_train)).^2);
    
    opt_corr_test(k) = corr(y_hat_test,y_test);
    opt_corr_train(k) = corr(y_hat_train,y_train);

end 
   
%--------------- End outer loop ---------------%

% ------------------------------------------------------------ 
% Prepare output data   
% ------------------------------------------------------------ 
    
% ADD OPTIMAL.YHAT

% Prediction performance 
% across the test sets 
optimal.bic_test =      opt_bic_test;
optimal.mse_test =      opt_mse_test;
optimal.nmse_test = 	opt_nmse_test;
optimal.corr_test = 	opt_corr_test;

% Prediction performance 
% across the train sets 
optimal.bic_train =     opt_bic_train;
optimal.mse_train =     opt_mse_train;
optimal.nmse_train =    opt_nmse_train;
optimal.corr_train =    opt_corr_train;

% Estimated model
% across folds 
optimal.efp =       opt_coef;
optimal.lambda =    opt_lambda;
optimal.rho =       opt_rho;
optimal.df =        opt_df;

end