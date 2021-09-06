function [ses_coef, T, optimal] = kfold_one_class_par(EEG, BOLD, path_pars, path_hc, varargin)

%   [tree, optimal] = kfold_cv_sessopns_par(EEG, BOLD,...) performs  
%   hierarchical clustering and nested k-fold cross-validation to obtain
%   a one-class model for the entire dataset
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
%     'rho'            The complexity paramter rho 
%     'lambda'         The complexity parameter lambda
%   
%   Outputs:
%
%     ses_coef         The coefficients of the model for each session 
%     T                The cluster assignment 
%     optimal          The struct with information corresponding  
%                      the best across each of the leave-one-out CV folds 
%       'rho'          The regularization parameter rho of the L21+1 fit
%       'lambda'       The complexity parameter lambda of the L21+1 fit 
%       'efp'          The activation pattern of the final fit
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

% Number of samples (time-points)
n_pnts = size(BOLD, 1);

% Number of regressors in the model
n_features = size(EEG, 2);

% Number of sessions 
n_sessions = size(EEG, 3);

% ------------------------------------------------------------ 
% Sanity check and process optional parameters 
% ------------------------------------------------------------ 

% Assign default values for each optional parameter
pnames = {'hc_distance' 'v' 'n' 'rho' 'lambda'}; 
dflts  = {["average", "correlation"] 0.3 10 [] []};

% Assign variables corresponding to optional parameters 
[hc_distance, V, N, rho, lambda] = internal.stats.parseArgs(pnames, dflts, varargin{:});

% Check if user supplied lambda and rho
user_supplied_rho = true;
user_supplied_lambda = true;

if isempty(rho)
    user_supplied_rho = false;
end   

if isempty(lambda)
    user_supplied_lambda = false;
end

if ~isequal(user_supplied_rho, user_supplied_lambda)
    error('Must supply both rho and lambda, or neither');
end

% ------------------------------------------------------------ 
% Hierarchical clustering 
% ------------------------------------------------------------ 

% First compute a model for each of the input sessions 

% In case the user hasn't supplied rho and lambda (must 
% supply both or neither), fix rho and compute the set 
% of 20 best lambda values for a random dataset
if ~user_supplied_rho   
    
    rho = 0.6;
    
    % Retreive pre-established lambda range of interest
    load(fullfile(path_pars, 'Lambdas.mat'), 'lambda');
    n_pars = length(lambda);
    
end

% EEG = reshape(EEG, [size(EEG, 1) size(EEG, 2) 3 23]);
% EEG = permute(EEG, [1 3 2 4]);
% EEG = reshape(EEG, [size(EEG, 1)*size(EEG, 2) size(EEG, 3) size(EEG, 4)]);
% 
% BOLD = reshape(BOLD, [size(BOLD, 1) 3 23]);
% BOLD = reshape(BOLD, [size(BOLD, 1)*size(BOLD, 2) size(BOLD, 3)]);
% 
% n_pnts = size(BOLD, 1);
% n_features = size(EEG, 2);
% n_sessions = size(EEG, 3);
 
if exist(fullfile(path_hc, 'ses_coef.mat'), 'file')
    
    load(fullfile(path_hc, 'ses_coef.mat'), 'ses_coef');
    
else
    
    % Pre-allocate variables 
    ses_coef = zeros(n_features + 1, n_sessions);

    % Go through all sessions and estimate a model for each 
    parfor se = 1 : n_sessions

        X = EEG(:, :, se);
        y = BOLD(:, se);

        % L2+1 elastic-net regression fit 
        [betas, stats] = lasso(X, y, 'Alpha', rho, 'Lambda', lambda);

        intercept = stats.Intercept;
        df = stats.DF';

        % Compute BIC for all lambda values and find minimum BIC 
        y_hat = intercept + X*betas; 
        mses = sum((y_hat - y).^2)'; 
        %nmse =  mses/sum((y - mean(y)).^2);  
        bic = log(length(y)).*df + length(y).*log(mses ./ length(y));
        [~, ind] = min(bic);
        ses_coef(:, se) = [intercept(ind); betas(:, ind)]; 

    end 
   
end

% Report clustering results 
report_clustering_results(ses_coef', hc_distance, fullfile('IMAGES', path_hc));
    
% Create a hierarchical cluster tree, Z
% Uses euclidean distance as the distance between leaves
% Uses the centroid of each cluster to compute distance 
% between clusters 
Z = linkage(ses_coef', hc_distance(1), hc_distance(2));

% Find the optimal number of clusters by finding the elbow 
% (or knee) method - find the knee of the objective 
% function, which in this case is the correlation function 
[n_clusters, ~] = knee_pt(Z(:, 3));

% Define clusters from the hierarchical cluster tree Z
T = cluster(Z, 'MaxClust', n_clusters);

% Find the largest cluster
[occurences, clusters]= hist(T, unique(T));
[~, ind] = max(occurences);
largest_cluster = clusters(ind);
 
% Select sessions belonging to the largest_cluster
sessions_in = (T == largest_cluster);
sessions_out = (T ~= largest_cluster);

% ------------------------------------------------------------ 
% Estimate one-class model and its performance 
% ------------------------------------------------------------ 

X = EEG(:, :, sessions_in);
y = BOLD(:, sessions_in);
K = size(X, 3);

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

%-------------- Begin outer loop ---------------%

% Leave-one-out CV procedure 
for k = 1 : K

    % Train set indices 
    ind_train = find(1 : K ~= k);
    
    % Test set indices 
    ind_test = k;
    
    % Train set
    X_train = X(:, :, ind_train); 
    X_train = reshape(permute(X_train, [1 3 2]), [n_pnts*(K-1) n_features]);
    y_train = y(:, ind_train); y_train = y_train(:);
    siz_train = length(y_train);
    
    % Test set
    X_test = squeeze(X(:, :, ind_test)); 
    y_test = squeeze(y(:, ind_test));
    siz_test = length(y_test);
    
    % Allocate bic and mse matrices for the learn and 
    % for the val set, each inner iteration through cols
    bic_val = zeros(n_pars, N); nmse_val = bic_val;
    df_inner = bic_val;

    %-------------- Begin inner loop --------------%  

    % The inner loop 
    % has N iterations 
    parfor n = 1 : N
        
        % Circular holdout assignment 
        ind_start = randi([1 siz_train]);
        indices_in = sort(crossvalind('Holdout', siz_train, V));
    
        % Assign learning set indices 
        ind_learn = zeros(size(indices_in));
        ind_learn(ind_start : end) = indices_in(1 : size(ind_learn) - ind_start + 1);
        ind_learn(1 : ind_start - 1) = indices_in(size(ind_learn) - ind_start + 2 : size(ind_learn));
        
        % Assign broadcast variables to loop variables for efficiency 
        % Large broadcast variables can cause significant communication 
        % between client and workers and increase parallel overhead
        % replace them for temporary variables, created inside the loop 
        X_par = X_train; y_par = y_train; lambda_par = lambda;
        
        % Assign validation set indices 
        ind_val = (~ind_learn);
        ind_val = find(ind_val);
        siz_val = length(ind_val);
        
        ind_learn = find(ind_learn); 
        
        % Assign learning and validation variables 
        X_learn = squeeze(X_par(ind_learn, :));
        y_learn = squeeze(y_par(ind_learn));
        X_val = squeeze(X_par(ind_val, :)); 
        y_val = squeeze(y_par(ind_val));

        % L2+1 elastic-net regression fit 
        [betas, stats] = lasso(X_learn, y_learn, ...
            'Alpha', rho, 'Lambda', lambda_par);
    
        intercept = stats.Intercept;
        df = stats.DF;
        
         df = flip(df)'; 
         betas = flip(betas, 2);
         intercept= flip(intercept);

        % Compute BIC for all lambda values and find minimum BIC 
        y_hat_val = intercept + X_val*betas; 
        mse_val = sum((y_hat_val - y_val).^2)'; 
        nmse_val(:, n) = mse_val/sum((y_val - mean(y_val)).^2);     
        bic_val(:, n) = log(siz_val).*df + siz_val.*log(mse_val ./ ...
            siz_val); 
        df_inner(:, n) = df;  
        
    end
    
    %--------------- End inner loop ---------------%    
    
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
    % all N inner iterations in the validation set
    aux_val(trash_lam)=Inf; [~, ind_opt] = min(aux_val); 
    if isempty(ind_opt); [~, ind_opt] = min(sum(bic_val, 2)); end
       
    % Save rho-lambda pair that minimizes
    % bic in learn and val sets combined 
    opt_lambda(k, 1) = lambda(ind_opt);
    opt_rho(k, 1) = rho;

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

% Session number 
optimal.ses =           sessions_in; 

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