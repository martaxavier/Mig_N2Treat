function [B,stats] = regress_L2_1_fast_greedy(X,Y,varargin)

%   [B,STATS] = regress_L2_1_fast_greedy(X,Y,...) 
%   performs L2 and L1-constrained linear least squares fits -
%   elastic net fit - relating the predictors in X to the responses 
%   in Y. In this version (i.e, "fast") the optimization method 
%   performs semi-greedy accelerated coordinate updates. 
%
%   Positional parameters:
%
%     X                A numeric matrix, NxP
%     Y                A numeric vector of length N
%   
%   Optional input parameters:  
%
%     'Alpha'          The ratio of L1 to L2 regularization incurred 
%     'NumLambda'      The number of lambda values to use, if the parameter
%                      'Lambda' is not supplied. Ignored if 'Lambda' is 
%                      supplied. regress_L2_1 may return fewer fits than
%                      specified by 'NumLambda' if the residual error of
%                      the fits drops below a threshold percentage of the 
%                      variance of Y.
%     'LambdaRatio'    Ratio between the minimum value and maximum value of
%                      lambda to generate, if the  parameter 'Lambda' is not 
%                      supplied. Legal range is [0,1[.  
%     'Lambda'         Complexity parameter of the L2+1 normalization.
%                      Will be returned in return argument STATS in 
%                      ascending order. The default is to generate a s
%                      sequence of lambda values, based on 'NumLambda'
%                      and 'LambdaRatio'. Regress_L2_1 will generate a 
%                      sequence, based on the values in X and Y, such that
%                      the largest lambda value is just sufficient to 
%                      produce all zero coefficients B. You may supply a 
%                      vector of real, non-negative values of lambda, in 
%                      place of its default sequence. 
%     'Standardize'    Whether to scale X prior to fitting the model
%                      sequence. This affects whether the regularization is
%                      applied to the coefficients on the standardized
%                      scale or the original scale. The results are always
%                      presented on the original data scale. Default is
%                      TRUE, do scale X.
%                      Note: X and Y are always centered.
%     'RelTolB'        Convergence threshold for coord descent algorithm
%                      The coordinate descent iterations will terminate
%                      when the relative change in the size of the
%                      estimated coefficients B drops below this threshold.
%     'RelTolMSE'      Convergence threshold for coord descent algorithm
%                      The coordinate descent iterations will terminate
%                      when the relative change in the size of the predi 
%                      ction MSE drops below this threshold.
%     'MaxIter'        Maximum number of iterations allowed.
%   
%   Return values:
%     B                The fitted coefficients for each model. 
%                      B will have dimension PxL, where 
%                      P = size(X,2) is the number of predictors, and
%                      L = length(lambda).
%     STATS            STATS is a struct that contains information about the
%                      sequence of model fits corresponding to the columns
%                      of B. STATS contains the following fields:
%
%       'Intercept'    The intercept term for each model. Dimension 1xL.
%       'Lambda'       The sequence of lambda penalties used, in ascending order. 
%                      Dimension 1xL.
%       'Alpha'        The alpha value that was used.
%       'DF'           The number of nonzero coefficients in B for each
%                      value of lambda. Dimension 1xL.
%       'MSE'          The mean squared error of the fitted model for each
%                      value of lambda. If cross-validation was performed,
%                      the values for 'MSE' represent Mean Prediction
%                      Squared Error for each value of lambda, as calculated 
%                      by cross-validation. Otherwise, 'MSE' is the mean
%                      sum of squared residuals obtained from the model
%                      with B and STATS.Intercept.
%

% --------------------------------------------------------------------
% Sanity check the positional parameters X,Y
% --------------------------------------------------------------------

% X a real 2D matrix
if ~ismatrix(X) || length(size(X)) ~= 2 || ~isreal(X)
    error('X is not a real 2D matrix');
end

if size(X,1) < 2
    error('Too few observations');
end

% Y a vector, same length as the columns of X
if ~isvector(Y) || ~isreal(Y) || size(X,1) ~= length(Y)
    error('Y is not a conforming vector');
end

% If Y is a row vector, convert 
% it to a column vector
if size(Y,1) == 1
    Y = Y';
end

% This screen (okrows) selects all the 
% predictions and response we can use
okrows = all(isfinite(X),2) & all(isfinite(Y),2);

% Remove observations with NaNs and 
% Infs in the predictor or response
X = X(okrows,:);
Y = Y(okrows);

% We need at least two observations after stripping NaNs and Infs.
if size(X,1) < 2
    error('Too few observations after removing NaN values');
end

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations. The corresponding coefficients
% will be returned as zero
constantPredictors = (range(X)==0);

% Index of non-constant predictors 
ever_active = ~constantPredictors;

% --------------------------------------------------------------------
% Parse and process the optional parameters
% --------------------------------------------------------------------

% Default value for 
% lambda ratio
lr_default = 2e-3;

% Default values for
% relative tolerances 
rtb_default = 1e-5; 
rtmse_default = 1e-5; 

% Default value for max
% number of iterations 
mi_default = 5e2; 

% Assign default values for each optional parameter
pnames = {'alpha' 'numlambda' 'lambdaratio' ...
    'lambda' 'standardize' 'reltolB' 'reltolMSE'  'maxiter'};
dflts  = { 1 100 lr_default [] true rtb_default rtmse_default mi_default};

[alpha, nLambda, lambdaRatio, lambda, ...
     standardize, reltolB, reltolMSE, maxIter] ...
     = internal.stats.parseArgs(pnames, dflts, varargin{:});

% === 'alpha' parameter ===

% Require 0 < alpha <= 1
if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || ...
        alpha < 0 || alpha > 1
    error('Invalid alpha')
end

% === 'Standardize' option ===

% Require a logical value
if ~isscalar(standardize) || (~islogical(standardize)...
        && standardize~=0 && standardize~=1)
    error('InvalidStandardize')
end

% === 'MaxIter' parameter ===

% mfilename is the name of the currently executing file
validateattributes(maxIter, {'numeric'},...
    {'scalar','positive','finite','integer'},...
    mfilename,'''MaxIter'' parameter');

% === 'Lambda' sequence ===

% lambdaMax is the penalty term (lambda) beyond which coefficients
% are guaranteed to be all zero.  If the command line does not provide
% a lambda sequence, we use lambdaMax in constructing the default 
% lambda sequence.  We always skip computation with lambda > lambdaMax
% because we know a priori that the computed coefficients will be zero
%
% nullMSE is the mse of the fit using just a constant term
% It is used to terminate the (ever-less penalized) fits when
% it becomes clear that we are overfitting

[lambdaMax, nullMSE] = computeLambdaMax(X, Y, alpha, standardize);

% Used with nullMSE (calculated below) to terminate
% (ever-less penalized) fits when overfitting is detected
userSuppliedLambda = true;

if isempty(lambda)
    
    % Used with nullMSE (calculated below) to terminate 
    % (ever-less penalized) fits when overfitting is detected
    userSuppliedLambda = false;
    
    % Sanity-check of 'NumLambda', should be positive integer
    if ~isreal(nLambda) || ~isfinite(nLambda) || nLambda < 1
        error('InvalidNumLambda');
    else
        nLambda = floor(nLambda);
    end
    
    % Sanity-checking of LambdaRatio, should be in [0,1[
    if ~isreal(lambdaRatio) || lambdaRatio <0 || lambdaRatio >= 1
        error('InvalidLambdaRatio');
    end
    
    if nLambda==1
        lambda = lambdaMax;
    else
        % Fill in a number "nLambda" of 
        % smaller values, on a log scale
        if lambdaRatio==0
                lambdaRatio = lr_default;
                addZeroLambda = true;
        else
            addZeroLambda = false;
        end
        lambdaMin = lambdaMax * lambdaRatio;
        loghi = log(lambdaMax);
        loglo = log(lambdaMin);
        lambda = exp(linspace(loghi,loglo,nLambda));
        if addZeroLambda
            lambda(end) = 0;
        else
            lambda(end) = lambdaMin;
        end
    end
    
else

    % Sanity check on user-supplied lambda
    % sequence. Should be non-neg real.
    if ~isreal(lambda) || any(lambda < 0)
        error('NegativeLambda');
    end

    lambda = sort(lambda(:),1,'descend');
    
end

% === 'RelTol' parameter ===

% Require 0 < RelTol < 1
if ~isscalar(reltolB) || ~isreal(reltolB) || ...
        ~isfinite(reltolB) || reltolB <= 0 || reltolB >= 1
    error('InvalidRelTol');
end


% --------------------------------------------------------------------
% Model fits
% --------------------------------------------------------------------

% The struct 'stats' will comprise the second return argument
% Put place holders for ever-present fields to secure the order
% we want in the struct
stats = struct();
stats.Intercept      = [];
stats.Lambda         = [];
stats.alpha            = alpha;
stats.DF             = [];
stats.MSE            = [];
stats.MSEs           = [];
stats.reachedMaxIter = [];

[B,Intercept,lambda,mse,mseIter,reachedMaxIter] = ...
    regressFit(X,Y,lambda,alpha,standardize,reltolB,...
    reltolMSE,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter);

% Store the number of non-zero
% coefficients for each lambda
df = sum(B~=0,1);

% --------------------------------------------------------------------
% Order results by ascending lambda
% --------------------------------------------------------------------

nLambda = length(lambda);
reverseIndices = nLambda:-1:1;
lambda = lambda(reverseIndices);
lambda = reshape(lambda,1,nLambda);
B = B(:,reverseIndices);
Intercept = Intercept(reverseIndices);
df = df(reverseIndices);
mse = mse(reverseIndices);
mseIter = mseIter(:,reverseIndices);
reachedMaxIter = reachedMaxIter(reverseIndices);

stats.Intercept = Intercept;
stats.Lambda = lambda;
stats.DF = df;
stats.MSE = mse;
stats.MSEs = mseIter;
stats.reachedMaxIter = reachedMaxIter;
   
end % regress_L2_1

% --------------------------------------------------------------------
% SUBFUNCTIONS 
% --------------------------------------------------------------------

% ===================================================
%                 regressFit() 
% ===================================================

function [B,Intercept,lambda,mspe,mseIter,reachedMaxIter] = ...
    regressFit(X,Y,lambda,alpha,standardize,reltolB,...
    reltolMSE,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter)
%
% ------------------------------------------------------
% Perform model fit for each lambda and the given alpha
% ------------------------------------------------------

[~,P] = size(X);
nLambda = length(lambda);

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ever_active & ~constantPredictors;

% === Standardtize variables ===

if standardize
    % Center and scale X
    [X0,muX,sigmaX] = zscore(X,1);
    % Avoid divide by zero with constant predictors
    sigmaX(constantPredictors) = 1;
else
    % Center X
    muX = mean(X,1);
    X0 = bsxfun(@minus,X,muX);
    sigmaX = 1;
end

% Center Y
muY = mean(Y);
Y0 = bsxfun(@minus,Y,muY);

% === Define the step size ===

% Compute the Lipschitz constant of each 
% coordinate - squared columnwise 2-norm
% of the least squares objective function, f
L = vecnorm(X0,2).^2;

% Compute the convexity parameter mu of 
% the least squares objective function, f
[~,S,~] = svd(X'*X); mu = S(end);
%mu = 0.01/N;

% Compute the step size
% for each coordinate 
gamma = 2./L;

% === Pre-allocate matrices ===

% b is be the current coefficient
% estimate, iteratively updated
% Because b is retained from one lam
% to the next, we get a warm start
b = zeros(P,1);

% Preallocate the returned matrix of
% coefficients, B, and the intercepts
B = zeros(P,nLambda);

% Preallocate the returned flag recahedMaxIter
reachedMaxIter = false(1,nLambda);

% Preallocate matrix of mse through iterations
mseIter = zeros(maxIter,nLambda);

% 1-by-P matrix of logical zeros 
active = false(1,P);

% === Go through lambdas ===

for i = 1:nLambda
    
    % Define lambda 
    lam = lambda(i);
    if lam >= lambdaMax
        continue;
    end
    
    % Initiate a counter
    % for the iterations
    ii = 1; 
    
    % Pre-allocate matrices with mse 
    % and b values for all iterations
    mse_ii = zeros(maxIter,1);
    b_ii = zeros(P,maxIter);
    
    % Coordinate descent 
    % initialization 
    t = sqrt(mu)/(1+mu); 
    w = mu*t;
    z = b;
    
    % Initiate animated line 
    %an = animatedline;
            
    % Iterative coordinate descent until
    % converged or reaches max iterations 
    for numIter = 1:maxIter
        
        bold = b;

        % Perform P iterations of coordinate descent 
        [b,z,active] = cdescentCycle(X0,Y0,b,z,t,w,...
            alpha,lam,gamma,active);
        
        % Save the mean squared error of current iteration, 
        % with the corresponding vector of coefficients b 
        bsig_ii = b ./ sigmaX';
        fit_ii = [ones(size(X,1),1) X] * ...
            [(muY-muX*bsig_ii); bsig_ii];
        residuals_ii = bsxfun(@minus, Y, fit_ii);
        mse_ii(ii) =  mean(residuals_ii.^2); 
        b_ii(:,ii) = b; 
        
        % Draw animated line 
        %addpoints(an,ii,mse_ii(ii));
        %addpoints(an,ii,norm( (b-bold)...
        %    ./ (1.0 + abs(bold)), Inf ));
        %drawnow limitrate
         
        % Update counter  
        % of iterations 
        ii = ii + 1; 
        
        % Check for convergence, in terms of coefficients or mse 
        if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) < reltolB || ...
            (ii > 30 && std(mse_ii(ii-30:ii-1))/...
            sqrt(mean(mse_ii(ii-30:ii-1))) < reltolMSE)
        
            % Cycling over the active set converged
            % Do one full pass through the predictors
            % If there is no predictor added to the active
            % set, we're done. Otherwise, resume coordinate 
            % descent iterations. 
            bold = b;
            potentially_active = thresholdScreen(X0, Y0, ...
                b, active, gamma(1)*lam*alpha);
            
            if any(potentially_active)
                new_active = active | potentially_active;
                [b,z,new_active] = cdescentCycle(X0,Y0,...
                    b,z,t,w,alpha,lam,gamma,new_active);
            else
                new_active = active;
            end

            if isequal(new_active, active)
                mseIter(:,i) = mse_ii;
                break
            else
                active = new_active;
            end
        
            if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) < reltolB || ...
            (ii > 30 && std(mse_ii(ii-30:ii-1))/...
            sqrt(mean(mse_ii(ii-30:ii-1))) < reltolMSE)
        
                mseIter(:,i) = mse_ii;
                break
                
            end
            
        end
        
        if numIter == maxIter
            
            warning('Maximum number of iterations reached');
            
            % In case the model doesn't converge, save the vector of 
            % coefficients that yielded the lowest mean squared error 
            % and use that as the final vector of coefficients 
            [~,idx]=min(mse_ii);b = b_ii(:,idx);
            
            reachedMaxIter(1,i) = true; 
            
        end

    end
    
    B(:,i) = b;
    mseIter(:,i) = mse_ii;
    
    % Halt if we have exceeded a threshold on the
    % percent of residual variance left unexplained
    if ~userSuppliedLambda
        
        % Calculate mse of the current fit
        bsig = b ./ sigmaX';
        fit = [ones(size(X,1),1) X] * [(muY-muX*bsig); bsig];
        residuals = bsxfun(@minus, Y, fit);
        mspe = mean(residuals.^2);

        if mspe < 1.0e-3 * nullMSE
            lambda = lambda(1:i);
            mseIter = mseIter(:,1:i);
            B = B(:,1:i);
            break
        end
        
    end
    
end % of lambda sequence

% ------------------------------------------
% Unwind the centering and scaling (if any)
% ------------------------------------------

B = bsxfun(@rdivide, B, sigmaX');
B(~ever_active,:) = 0;
Intercept = muY-muX*B;

% ------------------------------------------
% Calculate Mean Prediction Squared Error
% ------------------------------------------

BwithI = [Intercept; B];
fits = [ones(size(X,1),1) X]*BwithI;
residuals = bsxfun(@minus, Y, fits);
mspe = mean(residuals.^2);

if ~exist('mseIter','var')
    mseIter = zeros(1,nLambda);
end

end %-lassoFit

% ===================================================
%                 cdescentCycle() 
% ===================================================

function [b,z,active] = cdescentCycle(X0,Y0,b,z,t,w,...
            alpha,lam,gamma,active)

%    [b,z,t,active] = cdescentCycle(...) solves min_b f(b) + g(b)
%    where f is a smooth convex proper function and g is a convex
%    proper function with an easy to compute proximal operator
%
%    INPUTS:
%
%       X0          the predictior matrix 
%       Y0          the response variable 
%       b           the initial vector of coefficients
%       z           the initial auxiliar vector of coefficients 
%       t           the initial relaxation parameter
%       alpha       the ratio of L1/L2 regularization
%       lam         the complexity parameter
%       gamma       the step size 
%       active      the initial active set 
%
%    OUTPUTS:
%
%       b           the final vector of coefficients 
%       z           the final auxiliar vector of coefficients 
%       active      the final active set 

[N,~]=size(X0);

% Compute residuals
% for the current b 
r = Y0 - X0*b;

% Choose a random coord
% from the active set 
a = find(active);
a = a(randperm(numel(a)));

% Define the coordinate-wise 
% proximity operator of g = L_2+1 
proxg = @(xj,gamma)proxL2_1...
    (xj, alpha, lam, gamma);

% Define partial derivative of f
gradf = @(j,r)(- sum(X0(:,j).*r) / N); 

% Notation - j is the subscript 
% of the coordinate being updated,
% the subscript new corresponds to 
% the iteration k+1, and no subscript
% corresponds to the iteration k 

for j2 = a
    

    y = (1 - t).*b + t.*z;
        
    % Define j1 for current iteration
    [~,j1] = max(abs(X0'*(X0*y-Y0)));
    
    gammaj1 = gamma(j1);
    gammaj2 = gamma(j2);
    
    yj1 = y(j1); yj2 = y(j2);
   
    r = r - X0(:,j1)*(y(j1)-yj1);     
    bj_new = proxg( yj1 - gammaj1*gradf(j1,r), gammaj1);   
    bnew = b; bnew(j1) = bj_new;
    
    r = r - X0(:,j2)*(y(j2)-yj2);  
    u = ( t^2/(t^2+w) ).*z + ( w/(t^2+w) ).*y; uj2 = u(j2); 
    zj_new = proxg( uj2 - gammaj2*gradf(j2,r), gammaj2);   
    znew = z; znew(j2) = zj_new; 
    
    if bnew(j1) == 0
        active(j1) = false;
    end
    
    % Update variables for next cycle
    z = znew; b = bnew;
    
end


end %-cdescentCycle

% ===================================================
%                 proxL21_1() 
% ===================================================
function proxj = proxL2_1(xj, alpha, lam, gamma)
% proximal operator for  lambda(alpha||.||1 + (1-alpha)||.||2)

shrink_factor = 1 + gamma*lam*(1-alpha);
threshold = gamma*lam*alpha;
proxj = sign(xj) .* max( abs(xj) - threshold,0) .* shrink_factor;

end %proxL2_1

% ===================================================
%                 thresholdScreen() 
% ===================================================

function potentially_active = thresholdScreen(X0, Y0, ...
    b, active, threshold)

r = Y0 - X0(:,active)*b(active);

% We don't need the (b.*wX2)' term that one might expect
% because it is zero is zero for the inactive predictors
potentially_active = abs(r' *X0) > threshold;

end %-thresholdScreen


% ===================================================
%                 computeLambdaMaX() 
% ===================================================

function [lambdaMax, nullMSE] = computeLambdaMax(X, Y, alpha, standardize)
%
% lambdaMax is the penalty term (lambda) beyond which coefficients
% are guaranteed to be all zero.
%
% nullMse is the mse of the fit using just a constant term.
% It is provided in this function as a convenience, because it needs 
% to be calculated in the same context as lambdaMax whenever
% lambdaMax is calculated.

[N,~] = size(X);

% If we were asked to standardize the predictors, do so here because
% the calculation of lambdaMax needs the predictors as we will use
% them to perform fits.

if standardize
    % Center and scale X
    [X0,~,~] = zscore(X,1);
else   
    % Center X
    muX = mean(X,1);
    X0 = bsxfun(@minus,X,muX);    
end

% Center Y
muY = mean(Y);
Y0 = bsxfun(@minus,Y,muY);

% Calculate max lambda that  
% allows non-zero coefficients
dotp = abs(X0' * Y0);
lambdaMax = max(dotp) / (N*alpha);
nullMSE = mean(Y0.^2);

end
