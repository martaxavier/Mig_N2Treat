function [bestpairs,stats] = choose_rho_lambda(X,Y,siz_X,n_best,method,conv)

% Chooses best n_best rho-lambda pairs for input dataset and method 
% 
%   INPUTS:
%
%     X                  numeric matrix, NxP
%     Y                  numeric vector of length N
%     siz_X              the size of the problem
%     n_best             the number of best pairs to return
%     method             the method used for regression
%     conv               logical, determines if convergence is necessary
%
%   OUTPUTS:
%
%       bestpairs        n_best best rho-lambda pairs 
%
%       stats            struct that contains information about the
%                        sequence of model fits corresponding to the
%                        the searched parameters. stats contains the
%                        following fields:
%
%       'rho'            rho values for each of the
%                        rho-lambda combinations tested  
%       'lambda'         lambda values for each of the
%                        rho-lambda combinations tested  
%       'mse'            mean squared error of the fits for 
%                        each of the combinations tested 
%       'nmse'           normalized mean squared error of the 
%                        fits for each of the combinations tested 
%       'df'             degrees of freedom of the fits for 
%                        each of the combinations tested 
%       'bic'            bayesian information criterion of the 
%                        fits for each of the combinations tested 

tic 

% --------------------------------------------
% Input information 
% --------------------------------------------

% Default input data 
n_rho = 9;              % # of rhos to be tested 
n_lambda = 50;          % # of lambdas to be tested
lambda_ratio = 1e-4;    % ratio of lambda_min to lambda_max
max_iter = 6e2;         % # of maximum iterations of coodesc cycle

% If Y is a row vector, 
% convert to a column vector
if size(Y,1) == 1
    Y = Y';
end

[N,~] = size(X);

% Center and scale X
[X0,~,~] = zscore(X,1);

% Center Y
muY = mean(Y);
Y0 = bsxfun(@minus,Y,muY);

% Verify if method is one of the supported 
if ~strcmp(method,'l21_1') && ...
    ~strcmp(method,'elasticnet')
    error('Input method is not supported')
end

% --------------------------------------------
% Obtain array of rho and lambda values 
% --------------------------------------------

% Calculate array of rho values to be tested 
rho_min = 0.1; rho_max = 0.9;
rho = linspace(rho_min,rho_max,n_rho);

% Calculate max lambda that 
% permits non-zero coefficients
dotp = abs(X0' * Y0);
lambda_max = max(dotp) ./ (N*rho);

% Calculate minimum lambda using 
% the ratio between max and min lambda
lambda_min = lambda_max.*lambda_ratio;

% Obtain matrix of lambda values to be tested 
loghi = log(lambda_max);
loglo = log(lambda_min);
lambda = zeros(n_lambda,n_rho);

for ii = 1 : n_rho
%    [lambda(:,ii)] = get_lambdas(X,Y,rho(ii),...
%        'dfmin',10,'method',method,'numlambda',...
%        n_lambda,'lambdaratio',lambda_ratio);
    lambda(:,ii) = exp(linspace...
        (loghi(ii),loglo(ii),n_lambda));
    lambda(end,ii) = lambda_min(ii);
end

% --------------------------------------------
% Model fits for each rho-lambda pair
% --------------------------------------------

% Define fraction of total number 
% of combinations to be searched 
frac = 2;

% Generate random linear 
% indices the confound matrix
siz = n_lambda*n_rho;
idxs = randperm(siz,floor(siz/frac));

% Obtain indices in the confound matrix 
[lambda_idxs,rho_idxs] = ind2sub([n_lambda,n_rho],idxs);

% Allocate mse, df, nmse and bic matrices 
nmse = zeros(length(rho_idxs),1); df = nmse; 
Y_hat = zeros(N,length(rho_idxs));
MSEs = zeros(max_iter,length(rho_idxs));
reachedmaxiter = zeros(length(rho_idxs),1);

% Go through rho-lambda pairs 
parfor ii = 1 : length(rho_idxs)
    
    % Create temporary variables to 
    % improve parallel pool efficiency
    rho_par = rho; lambda_par = lambda;
    
    % Read current rho-lambda pair
    r = rho_par(rho_idxs(ii));
    lam = lambda_par(lambda_idxs(ii),rho_idxs(ii));
    
    % Screen the input method
    % Only two options possible 
    if strcmp(method,'l21_1')
    
        % Fit the model with current rho-lambda
        [betas,fit] = regress_L21_1(X,Y,siz_X,...
            'Rho',r,'Lambda',lam,'MaxIter',max_iter);
    
    else
        
        % Fit the model with current rho-lambda
        [betas,fit] = lasso(X,Y,'Alpha',r,...
            'Lambda',lam,'MaxIter',max_iter);
        
    end
    
    % Save mse and df of the fit 
    nmse(ii) = fit.MSE;
    df(ii) = fit.DF;
    Y_hat(:,ii) = fit.Intercept + ...
        X*betas;
    
    % Save variables returned only by l21_1
    if strcmp(method,'l21_1')
        num_iter = length(fit.MSEs);
        MSEs_par = [fit.MSEs; zeros(max_iter-num_iter,1)];
        MSEs(:,ii) = MSEs_par;
        reachedmaxiter(ii) = fit.reachedMaxIter;
    end
    
end

    % Compute mse and bic values 
    mse = sum((Y_hat-Y).^2)';
    bic = log(N).*df + N.*log(mse./N);

% --------------------------------------------
% Choose best pairs   
% --------------------------------------------

% Sort by ascendent mse 
[~,sortedindices_prior] = sort(bic);

if conv
    idxs_conv = reachedmaxiter(sortedindices_prior)==0;
    sortedindices = sortedindices_prior .* idxs_conv;
    sortedindices = sortedindices(sortedindices>0);
else
    sortedindices = sortedindices_prior;
end

% Assign rho and lambda of best rho-lambda pairs 
bestpairs = zeros(n_best,2);
bestpairs(:,1) = rho(rho_idxs(sortedindices(1:n_best)));
bestpairs(:,2) = lambda(lambda_idxs(sortedindices(1:n_best)));

% --------------------------------------------
% Organize outputs  
% --------------------------------------------
    
    % Organize stats struc
    stats.rho = rho(rho_idxs(sortedindices_prior))';
    stats.lambda = lambda(lambda_idxs(sortedindices_prior))';
    stats.rho = stats.rho(sortedindices_prior);
    stats.lambda = stats.lambda(sortedindices_prior);
    stats.mse = mse(sortedindices_prior);
    stats.nmse = nmse(sortedindices_prior);
    stats.df = df(sortedindices_prior);
    stats.bic = bic(sortedindices_prior);   
    stats.reachedmaxiter = reachedmaxiter(sortedindices_prior);
    stats.summary = struct2table(stats);
  
 toc
 
end