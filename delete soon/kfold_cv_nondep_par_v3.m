function [model,optimal] = kfold_cv_nondep_par_v3(EEG,BOLD,varargin)

%   [model,optimal] = kfold_cv(EEG,BOLD,...) performs non-dependent,  
%   nested k-fold cross-validation to obtain the best model for the
%   input EEG and BOLD data
%   Version v3: pararelled computing; inner loop is also non-dependent
%
%   Input data:
%
%     EEG              The model's X 
%     BOLD             The model's Y
%   
%   Optional input parameters:  
%
%     'k'              The number of folds in the outer cv loop 
%     'v'              The number of folds in the inner cv loop 
%     'inneriter'      The number of iterations in the inner cv loop 
%     'regress'        The regression method used to fit the model
%                      'elasticnet' or 'l21_1' 
%     'autocorr'       The order of the auto-correlation function 
%     'rho'            Thecomplexity parameter rho of the L21+1 fit 
%     'lambda'         The complexity parameter lambda of the L21+1 fit
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

% X a real 2D matrix
if ~ismatrix(EEG) || length(size(EEG)) ~= 2 || ~isreal(EEG)
    error('EEG is not a real 2D matrix');
end

if size(EEG,1) < 2
    error('Too few observations');
end

% If Y is a row vector, convert to a column vector
if size(BOLD,1) == 1
    BOLD = BOLD';
end

% Y a vector, same length as the columns of X
if ~isvector(BOLD) || ~isreal(BOLD) || size(BOLD,1) ~= length(BOLD)
    error('BOLD is not a conforming vector');
end

% Number of samples (time-points)
n_pnts = length(BOLD);

% Number of regressors in the model
n_features = size(EEG, 2);


% ------------------------------------------------------------ 
% Sanity check and process optional parameters 
% ------------------------------------------------------------ 

% Assign default values for each optional parameter
pnames = {'k' 'v' 'inneriter' 'regress' 'rho' 'lambda',...
    'numpars' 'sizx' 'autocorr'}; dflts  = { 15 10 ...
    10 'l21_1' [] [] 20 [n_pnts 31*6 4] 2};

% Assign variables corresponding to optional parameters 
[K, V, N, method, rho, lambda, n_pars, siz_X, h] = ...
    internal.stats.parseArgs(pnames, dflts, varargin{:});

% Check if the method provided is 
% one of the methods supported
if ~strcmp(method,'l21_1') && ...
    ~strcmp(method,'elasticnet')
    error('Input method is not supported')
end

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

% Allocate arrays of optimal parameters
% for each test/train pair
opt_lambda = zeros(K,1); opt_rho = opt_lambda;
opt_coef = zeros(n_features+1,K);
opt_df = opt_lambda; opt_mse = opt_lambda;
opt_nmse = opt_lambda; opt_bic = opt_lambda; 
opt_corr = opt_lambda; siz_dep = opt_lambda;
opt_pval = opt_lambda;

% In case the user hasn't supplied rho and lambda, 
% (must supply both or neither), fix rho and compute
% the set of 20 best lambda values for a random data  
% set of the size of the train set. Use the computed
% values for all the train/test sets. This is done
% because the optimal rho-lambda depend largely on 
% the size of the data
if ~user_supplied_rho   
    
    rho = 0.6;
    [lambda] = get_lambdas(EEG,BOLD,siz_X,rho,...
        'dfmin',45,'method',method,...
        'numlambda',n_pars,'lambdaratio',1e-1); 
    
end

%-------------- Begin outer loop ---------------%

% Divide time-series into K total folds 
% (assign each time-point to one fold)
indices_out = crossvalind('Kfold',n_pnts,K);

for k = 1 : K

    % Assign test set indices 
    idx_test = (indices_out == k); 
    idx_test = find(idx_test);
    siz_test = length(idx_test);
    
    % Find samples correlated with all samples within test set
    % and remove dependent samples from the training set 
    rem = -h : h; dep = idx_test - rem; 
    dep = reshape(dep,[size(dep,1)*size(dep,2),1]);
    dep = unique(dep); dep(dep>n_pnts) = []; dep(dep<1) = [];
    siz_dep(k) = length(dep) - siz_test;
    
    % Assign training set indices 
    idx_train = (1:n_pnts)'; 
    idx_train(dep) = []; 
    siz_train = length(idx_train);
    siz_test = length(idx_test);
    
    % Assign train and test X (EEG) and Y (BOLD) variables 
    X_train = EEG(idx_train,:); y_train = BOLD(idx_train);
    X_test = EEG(idx_test,:); y_test = BOLD(idx_test);
    
    % Allocate bic and mse matrices for the learn and 
    % for the val set, each inner iteration through cols
    %bic_learn = zeros(n_pars,N); nmse_learn = bic_learn; 
    bic_val = zeros(n_pars,N); nmse_val = bic_val;
    df_inner = zeros(n_pars,N);

    %-------------- Begin inner loop --------------%

    % The inner loop 
    % has N iterations 
    parfor n = 1 : N
    
        % Assign broadcast variables to loop variables for efficiency 
        % Large broadcast variables can cause significant communication 
        % between client and workers and increase parallel overhead
        % replace them for temporary variables, created inside the loop 
        EEG_par = EEG; BOLD_par = BOLD; lambda_par = lambda;

        % Divide train set into K equally sized folds 
        indices_in = crossvalind('Kfold',siz_train,V);
    
        % Assign test set indices 
        idx_val = (indices_in == 1); 
        idx_val = find(idx_val);
        siz_val = length(idx_val);
        
        % Remove dependencies 
        rem = -h : h; dep = idx_val - rem; 
        dep = reshape(dep,[size(dep,1)*size(dep,2),1]);
        dep = unique(dep); dep(dep>siz_train) = []; dep(dep<1) = [];
        
        % Assign learning set indices 
        idx_learn = (1:siz_train)'; 
        idx_learn(dep) = []; 
        %siz_learn = length(idx_learn);
        
        % Assign learning and validation variables 
        X_learn = EEG_par(idx_learn,:);y_learn = BOLD_par(idx_learn);
        X_val = EEG_par(idx_val,:);y_val = BOLD_par(idx_val);
            
        % Screen the input method
        if strcmp(method,'l21_1')

           % L21+1 fit with rho = r and lambda = l
            [betas,stats] = regress_L21_1(X_learn,...
            y_learn,siz_X,'Rho',rho,'Lambda',lambda_par,'MaxIter',1e3);
            [~,col] = find(betas); df = accumarray(col,1); 
            df(setdiff(1:n_pars,col))=0; df = flip(df);
            betas = flip(betas,2);

        else

            % L2+1 fit with rho = r and lambda = l
            [betas,stats] = regress_L2_1(X_learn,...
            y_learn,'Alpha',rho,'Lambda',lambda_par,'MaxIter',1e3);
            df = stats.DF; df = flip(df)'; betas = flip(betas,2);

        end
        
        % Save intercept values for current rho-lambda
        intercept = stats.Intercept; intercept=flip(intercept);

        % Compute bic values for all rho-lambda 
        % pairs in the val and learn set 
        y_hat_val = intercept + X_val*betas;      
        mse_val = sum((y_hat_val-y_val).^2)';     
        nmse_val(:,n) = mse_val/sum((y_val - mean(y_val)).^2);     
        bic_val(:,n) = length(siz_val).*df + ...
            length(siz_val).*log(mse_val./siz_val);
        
        df_inner(:,n) = df;
        
    end
    
    %--------------- End inner loop ---------------%
    
    % Flag as ineligible, in 'trash_lam' (1) lambda values for which 
    % nmse values are, at any of the inner iterations, above a given 
    % inacceptable threshold, 'thresh_nmse' and (2) lambda values for 
    % which dof values are, at any of the inner iterations, zero
    % Flag, in 'trash_n', inner iterations in which the average nmse 
    % value in the validation set was above 'thresh_nmse' 
    thresh_nmse = 0.95; aux = nmse_val>thresh_nmse; 
    aux(:,mean(nmse_val,1)>thresh_nmse-0.05)=0; 
    trash_lam = find(sum(aux,2)>=1); 
    [rows,~] = find(df_inner==0); 
    trash_lam = unique([trash_lam; unique(rows)]);
    
    % Average bic values without considering
    % iterations flagged as "trash" 
    %trash_n = find(mean(nmse_val,1)>thresh_nmse);
    %aux = bic_val; aux(:,trash_n)=0; 
    aux = bic_val; aux(:,mean(nmse_val,1)>thresh_nmse)=0; 
    [rows,~,val]=find(aux);
    aux_val = accumarray(rows,val,[],@mean);
    
    % Find optimal rho-lambda for the current test/train pair,
    % i.e, find rho-lambda that minimizes sum of bics through  
    % all N inner iterations in the validation set
    aux_val(trash_lam)=Inf; [~, idx_opt] = min(aux_val); 
    if isempty(idx_opt); [~,idx_opt] = min(sum(bic_val,2)); end
       
    % Save rho-lambda pair that minimizes
    % bic in learn and val sets combined 
    opt_lambda(k,1) = lambda(idx_opt);
    opt_rho(k,1) = rho;

    if strcmp(method,'l21_1')

        % Compute L21+1 coefficients for the current iteration of  
        % the outer CV procedure (k), using optimal lambda parameter 
        [betas,stats] = regress_L21_1(X_train,y_train,siz_X,...
        'Rho',opt_rho(k),'Lambda',opt_lambda(k),'MaxIter',1e3);
        betas(abs(betas)<5e-4)=0; 
        opt_df(k)=length(find(betas));

    else

        % Compute L2+1 coefficients for the current iteration of  
        % the outer CV procedure (k), using optimal lambda parameter 
        [betas,stats] = regress_L2_1(X_train,y_train,...
        'Alpha',opt_rho(k),'Lambda',opt_lambda(k),'MaxIter',1e3);
        opt_df(k) = stats.DF; 

    end
    
    % Compute the y hat for the 
    % current outer iteration k
    y_hat = stats.Intercept + X_test*betas;
    
    % Compute the remaining parameters
    % for the current outer iteration k
    opt_coef(:,k) = [stats.Intercept; betas];
    opt_mse(k) = sum((y_hat-y_test).^2);
    
    opt_bic(k) = log(siz_test).*opt_df(k) + ...
        length(siz_test).*log(opt_mse(k)./siz_test);

    opt_nmse(k) = opt_mse(k)/sum((y_test - mean(y_test)).^2);
    [opt_corr(k),opt_pval(k)] = corr(y_hat,y_test);

end 
   
%--------------- End outer loop ---------------%

% ------------------------------------------------------------ 
% Estimate a final model  
% ------------------------------------------------------------ 

% Use the mode of the rho-lambda
% pairs across the K cv folds 
lam = mode(opt_lambda);
r = mode(opt_rho);

% Estimate the models final  
% set of coefficients 
if strcmp(method,'l21_1')
    
    [betas,stats] = regress_L21_1(EEG,...
    BOLD,siz_X,'Rho',r,'Lambda',lam,'MaxIter',1e3);
    betas(abs(betas)<5e-4)=0; 

else
    
    [betas,stats] = regress_L2_1(EEG,...
    BOLD,'Alpha',r,'Lambda',lam,'MaxIter',1e3);

end

model.efp = [stats.Intercept;betas];


% ------------------------------------------------------------ 
% Prepare output data   
% ------------------------------------------------------------ 

% Average model performance
% parameters across all K folds
% and assign model structure 
model.bic = mean(opt_bic);
model.mse = mean(opt_mse);
model.nmse = mean(opt_nmse);
model.corr = mean(opt_corr);

tval = model.corr*sqrt((siz_test-2)/(1-model.corr^2));
pval = 2*tcdf(tval,siz_test-2,'upper');
model.pval = pval;

model.lambda = lam;
model.rho = rho;
model.time = toc;

% Assign model parameters that depend
% on the final set of coefficients
model.df = length(find(model.efp));
model.yhat = model.efp(1) ...
    + EEG*model.efp(2:end);

% Assign optimal structure
optimal.bic = opt_bic;
optimal.mse = opt_mse;
optimal.nmse = opt_nmse;
optimal.corr = opt_corr;
optimal.efp = opt_coef;
optimal.df = opt_df;
optimal.lambda = opt_lambda;
optimal.rho = opt_rho;
optimal.pval = opt_pval;

end