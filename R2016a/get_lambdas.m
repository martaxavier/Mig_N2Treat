function [lambda] = get_lambdas(X,Y,siz_X,rho,varargin)

% Gets list of lambda ranges for input rho value and dataset
% 
%   INPUTS:
%
%     X                  numeric matrix, NxP
%     Y                  numeric vector of length N
%     rho                rho parameter value or array
%                        of rho parameter values 
%
%   Optional input parameters
%    
%     method             the method used for regression
%                        may be 'elasticnet' or 'l21_1'
%     dfmin              the minimum number of dof in the model
%     lambdaratio        the ratio between the minimum value and
%                        maximum value of lambda
%     numlambda          The number of lambda values to use
%
%   OUTPUTS:
%
%       lambda           matrix containing the array of 
%

% --------------------------------------------
% Check and process input parameters 
% --------------------------------------------

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

% If rho is a column vector,
% convert to a row vector 
if size(rho,2) == 1
    rho = rho';
end

n_rho = length(rho);

% --------------------------------------------
% Check and process optional parameters 
% --------------------------------------------

% Assign default values for each optional parameter
pnames = {'method' 'dfmin' 'lambdaratio' 'numlambda'};
dflts  = { 'elasticnet' 1 1e-3 20};

[method, df_min, lambda_ratio, n_lambda] ...
     = internal.stats.parseArgs(pnames, dflts, varargin{:});
 
 % Verify if method is one of the supported 
if ~strcmp(method,'l21_1') && ...
    ~strcmp(method,'elasticnet')
    error('Input method is not supported')
end

% --------------------------------------------
% Obtain array of lambda values 
% --------------------------------------------

% Calculate max lambda that 
% permits non-zero coefficients
dotp = abs(X0' * Y0);
lambda_max = max(dotp) ./ (N*rho);

df = 1; 
lambda = lambda_max - ...
    log(linspace(exp(0.1*lambda_max),exp(lambda_max),30));
lambda(end) = 0.01;
ii = 1;
lam = lambda(ii);

while df < df_min
    lam = lambda(ii);
    if strcmp(method,'l21_1')
       [~,stats]=regress_L21_1(X,Y,siz_X,'Rho',rho,'Lambda',lam,'MaxIter',6e2);
    else
       [~,stats]=regress_L2_1(X,Y,'Alpha',rho,'Lambda',lam,'MaxIter',6e2);
    end
    df = stats.DF;
    ii = ii + 1;
end 

lambda_max = lam;

% Calculate minimum lambda using 
% the ratio between max and min lambda
lambda_min = lambda_max.*lambda_ratio;

% Obtain matrix of lambda values to be tested 
loghi = log(lambda_max);
loglo = log(lambda_min);
lambda = zeros(n_lambda,n_rho);

for ii = 1 : n_rho
    lambda(:,ii) = exp(linspace...
        (loghi(ii),loglo(ii),n_lambda));
    lambda(end,ii) = lambda_min(ii);
end

end