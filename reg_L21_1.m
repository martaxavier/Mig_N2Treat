function [B, stats] = reg_L21_1(X, Y, siz_X, varargin)

% --------------------------------------------------------------------
% Read problem 
% --------------------------------------------------------------------

lam_ratio_dflt = 1e-3;
max_iter_dflt = 1e3;
tol_b_dflt = 1e-12;

var_names = {'alpha','numLam' 'lamRatio' 'lam' 'maxIter' 'tol'};
dflts = { 1 100 lam_ratio_dflt [] max_iter_dflt tol_b_dflt};

[alpha, n_lam, lam_ratio, lam, max_iter, tol_b] = ...
    internal.stats.parseArgs(var_names, dflts, varargin{:});

[lam_max, ~] = compute_lam_max(X, Y, alpha);
lam = get_lam(lam, n_lam, lam_max, lam_ratio, lam_ratio_dflt);

% --------------------------------------------------------------------
% Model fits
% --------------------------------------------------------------------

stats = struct();
stats.intercept         = [];
stats.lam            = [];
stats.alpha             = alpha;
stats.DF                = [];
stats.MSE               = [];
stats.MSEs              = [];
stats.reached_max_iter  = [];

[B, intercept, lam, mse, mse_iter, reached_max_iter] = ...
    regress_fit(X, Y, siz_X, lam, alpha, lam_max, max_iter, tol_b);

% Store the number of non-zero
% coefficients for each lam
df = sum(B~=0, 1);

% --------------------------------------------------------------------
% Order results by ascending lambda
% --------------------------------------------------------------------

n_lam = length(lam);
reverse_indices = n_lam : -1 : 1;
lam = lam(reverse_indices);
lam = reshape(lam, 1, n_lam);
B = B(:, reverse_indices);
intercept = intercept(reverse_indices);
df = df(reverse_indices);
mse = mse(reverse_indices);
mse_iter = mse_iter(:, reverse_indices);
reached_max_iter = reached_max_iter(reverse_indices);

stats.intercept = intercept;
stats.lam = lam;
stats.DF = df;
stats.MSE = mse;
stats.MSEs = mse_iter;
stats.reached_max_iter = reached_max_iter;
   
end % reg_L2_1

% --------------------------------------------------------------------
% SUBFUNCTIONS 
% --------------------------------------------------------------------

% ===================================================
%                 regress_fit() 
% ===================================================
function [B, intercept, Lam, mse, mse_all, reached_max_iter] = ...
    regress_fit(X, Y, siz_X, Lam, alpha, lam_max, max_iter, tol_b)
%
% ------------------------------------------------------
% Perform model fit for each lambda and the given rho
% ------------------------------------------------------

[~, p] = size(X);
n_lam = length(Lam);

% Center and scale X
[X0, muX, sigmaX] = zscore(X,1);

% Center Y
muY = mean(Y);
Y0 = bsxfun(@minus, Y, muY);

% Pre-allocate 
b = zeros(p, 1);
B = zeros(p, n_lam); 

%reached_max_iter = false(1, n_lam);
%mse_all = zeros(max_iter, n_lam);

for i = 1 : n_lam
    
    % Define lambda 
    lam = Lam(i);
    if lam >= lam_max
        continue;
    end
    
    X0 = reshape(X0, [size(X,1) siz_X]);
    b = reshape(b, siz_X);
    
    % Lipschitz constant
    if length(size(X0)) == 2
        L = norm(X0)^2;
    elseif length(size(X0)) == 3
        XX = [];
        for ii = 1:size(X0,1)
            XX = [XX squeeze(X0(ii,:,:))]; 
        end
        L = norm(XX)^2;
    end
    
    % Function F and G 
    F = @(x)lam*norm(x, 1);
    G = @(x)1/2*norm(Y0 - X0*x)^2;
    if length(size(X0)) == 3
        G = @(x) G3(X0, x, Y0);
    end

    % Function to record the energy
    options.report = @(x)F(x)+G(x);

    % Bench the algorithm
    options.niter = max_iter;

    % Specify method
    options.method = 'fista'; 

    % Proximal operator of F 
    ProxF = @(x, tau)prox_L21_1(x, lam*tau, lam*alpha*tau);

    % Gradient operator of G
    GradG = @(x)X0'*(X0*x - Y0);
    if length(size(X0)) == 3
        GradG = @(x) gradG3(X0, x, Y0);
    end

    % Forward-backward optimization 
    [b, ~] = perform_fb(b, ProxF, GradG, L, options);
    
    % Initiate a counter
    % for the iterations
%     iter = 1; 
    
    % Pre-allocate 
%     mse_iter = zeros(max_iter, 1);
%     b_iter = zeros(p, max_iter);
%     
    % Initiate animated line
    %an = animatedline;
            
    % Iterative coordinate descent until
    % converged or reaches max iterations 
%     for num_iter = 1 : max_iter
%         
%         b_old = b;
% 
%         % Save the mean squared error of current iteration, 
%         % with the corresponding vector of coefficients b 
%         b_sig_iter = b ./ sigmaX';
%         fit_iter = [ones(size(X,1),1) X] * [(muY-muX*b_sig_iter); b_sig_iter];
%         residuals_iter = bsxfun(@minus, Y, fit_iter);
%         mse_iter(iter) =  mean(residuals_iter.^2); 
%         b_iter(:, iter) = b; 
%         
%         % Draw animated line      
%         %addpoints(an,ii,norm( (b-bold) ...
%         %./ (1.0 + abs(bold)), Inf ));
%         %drawnow limitrate  
%         
%         % Update the counter
%         % for the iterations
%         iter = iter + 1; 
%         
%         % Check for convergence
%         if norm( (b - b_old) ./ (1.0 + abs(b_old)), Inf ) < tol_b
%             
%                 mse_all(:, i) = mse_iter;
%                 break
%                 
%         end
%         
%         if num_iter == max_ter
%             
%             warning('Maximum number of iterations reached');
% 
%             [~, idx] = min(mse_iter); b = b_iter(:, idx);
%             reached_max_iter(1, i) = true; 
%             
%         end
% 
%     end
    
%    mse_all(:, i) = mse_iter;
    B(:, i) = b(:);
    
end % of lambda sequence

% ------------------------------------------
% Unwind the centering and scaling (if any)
% ------------------------------------------

B = bsxfun(@rdivide, B, sigmaX');
intercept = muY - muX*B;

% ------------------------------------------
% Calculate Mean Prediction Squared Error
% ------------------------------------------

BwithI = [intercept; B];
fits = [ones(size(X,1),1) X]*BwithI;
residuals = bsxfun(@minus, Y, fits);
mse = mean(residuals.^2);

end % regress_fit 


% ===================================================
%                 prox_L21_1() 
% ===================================================

function [x,R] = perform_fb(x, ProxF, GradG, L, options)
%
%   Solves min_x g(x) + f(x)
%   where g is a smooth convex proper function and f is a
%   convex proper function with an easy to compute proximal operator
%
%   Use several first order-scheme depending on options.method:
%       options.method = 'fb' : classical Foward-backward
%       options.method = 'fista' : FISTA method of Beck and Teboule
%       options.method = 'nesterov' : Nesterov scheme
%
%   INPUTS:
%   ProxF(y,sigma) computes Prox_{sigma*F}(x)
%   GradG(x) computes \nabla f(x)
%   L is the lipschitz constant of the gradient, if g is C^2:
%       L = max_x norm( Hg(x) ) 
%       where Hg(x) is the hessian of g at x. 
%       For instance, if g(x)=1/2*|A*x-y|^2 then tau = norm(A)^2.
%   options.niter is the number of iterations.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

options.null = 0;
method = options.method;
report = options.report;
niter = options.niter;
verb = 0;
fbdamping = 1.8;

t = 1;  % fista & nesterov
tt = 2/L; gg = 0; A = 0; % nesterov
y = x;
x0 = x;
for i=1:niter 
  	R(i) = report(x);
    if verb
        progressbar(i,niter);
    end
    switch method
        case 'fb'
            x = ProxF( x-fbdamping/L*GradG(x), fbdamping/L );
        case 'fista'
            xnew = ProxF( y - 1/L*GradG(y), 1/L );
            tnew = (1+sqrt(1+4*t^2))/2;
            y = xnew + (t-1)/(tnew)*(xnew-x);
            x = xnew; t = tnew;
        case 'nesterov'
            a = (tt + sqrt(tt^2 + 4*tt*A))/2;
            v = ProxF( x0-gg, A );
            z = (A*x+a*v)/(A+a);
            x = ProxF( z - 1/L*GradG(z) , 1/L  );
            gg = gg +  a * GradG(x); % P'*(P*x-y);
            A = A + a;
        otherwise
            error('Unknown method');
            
    end      
end

end

% ===================================================
%                 prox_L21_1() 
% ===================================================

function prox = prox_L21_1(x, lambda, alpha)
 % proximal operator for  alpha||.||1 + lambda||.||21
 % x: ExF
 % mu: large value will lead to spatially very sparse solution
 % lambda: large value will promote sources with smooth time series
 % Proximal operator added by Claire Cury
 
 for p=1:size(x,1)
     for k = 1:size(x,2)
         div_=0;
         for kk = 1:size(x,2)
             div_ = div_ +max(abs(x(p,kk))-alpha,0)^2;
         end
         div_ = sqrt(div_);
         if div_ ==0
             frac =0;
         else
             frac = lambda/div_;
         end
         prox(p,k) = sign(x(p,k)) * max(abs(x(p,k))-alpha,0) * max((1 - frac), 0);
     end
 end
 
end
 
% ===================================================
%                 gradG3() 
% ===================================================

function gradG = gradG3(A, x, y)
% A: TxExF
% x: ExF
% y: Tx1
    gradG=zeros(size(x));
    for i=1:size(A,1)
        gradG=gradG+squeeze(A(i,:,:))*(trace(squeeze(A(i,:,:))'*x) - y(i));
    end
end

% ===================================================
%                 G3()
% ===================================================

function G = G3(A, x, y)
    G=0;
    for i=1:size(A,1)
        G=G+(y(i) - trace(squeeze(A(i,:,:))'*x))^2;
    end
end
        
% ===================================================
%                 compute_lam_max() 
% ===================================================

function [lam_max, null_mse] = compute_lam_max(X, Y, alpha)
%
% lam_max is the penalty term (lam) beyond
% which coefficients are guaranteed to be all zero

[N,~] = size(X);

% Center and scale X
[X0,~,~] = zscore(X,1);
    
% Center Y
muY = mean(Y);
Y0 = bsxfun(@minus,Y,muY);

% Calculate max lam that
% allows non-zero coefficients
dotp = abs(X0' * Y0);
lam_max = max(dotp) / (N*alpha);
null_mse = mean(Y0.^2);

end

% ===================================================
%                 get_lam() 
% ===================================================

function [lam] = get_lam(lam, n_lam, lam_max, lam_ratio, lam_ratio_dflt)
    
if isempty(lam)

    if n_lam == 1
        
        lam = lam_max;  
        
    else  
        
        if lam_ratio == 0
            lam_ratio = lam_ratio_dflt;
            add_zero_lam = true;  
        else
            add_zero_lam = false;
        end
    
        lam_min = lam_max * lam_ratio;
        loghi = log(lam_max);
        loglo = log(lam_min);
        lam = exp(linspace(loghi, loglo, n_lam));
        
        if add_zero_lam
            lam(end) = 0;
        else
            lam(end) = lam_min;
        end
        
    end
    
else
    
lam = sort(lam(:), 1, 'descend');

end

end
