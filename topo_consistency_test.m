function [decision, pval, but_one_prob] = topo_consistency_test(data, n_rand, metric, path_pars)
%
%   [prob] = topo_consistency_test tests, for the current metric, the 
%            consistency of the group-level topographies, through  
%            a global randomization test 
%                 
%   Inputs:
%
%     data          	    the topographic data
%     n_rand                the number of randomizations 
%     metric                the current metric 
% 
%    Outputs:
%
%    prob                   the probability of rejecting the 
%                           null-hypothesis 
%    prob_but_one           the probability of rejecting the 
%                           null-hypothesis, after the removal
%                           of each observation 
%

% Get parameters of current metric 
get_metric_pars;

% First-level stats 
data = data';

% Number of observations
n_obs = size(data, 1);

% Reshape data to be obs x channels x features 
data = reshape(data, [n_obs dim(1) dim(2)*dim(3)]);

% Test topographic consistency of the data 
[decision, pval] = global_rand_test(data, n_rand);

% Test the effect of removing one observation in 
% the topographic consistency of the data 
but_one_prob = zeros([n_obs size(decision)]);

% Test the impact that removing each
% of the observations has on the 
% estimated probability of significance 
% for o = 1 : n_obs
%    obs = 1 : n_obs;
%    obs(o) = [];
%    but_one_rho = ...
%        squeeze(data(obs, :, :, :, :));
%    but_one_prob(o, :, :) = ...
%        global_rand_test(but_one_rho, n_rand, dim);
% end

end

% ============================================================
% [decision, pval] = global_rand_test(data, n_rand)         
% ============================================================

function [decision, pval] = global_rand_test(data, n_rand)

    % Dimensions of input data 
    n_obs = size(data, 1);
    n_chans = size(data, 2);
    n_features = size(data, 3); 

    % Allocate and initialize global field power 
    gfp = zeros(n_rand + 1, n_features);
    gfp(1, :) = std(squeeze(mean(data, 1)), 1);
    
    % Allocate random data variable 
    data_rand = zeros(size(data));
    
    for r = 1 : n_rand
        
        for o = 1 : n_obs
            
            for f = 1 : n_features 
            
                % Compute random data for current observation by 
                % randomnizing the channel data 
                data_rand(o, :, f) = data(o, randperm(n_chans), f);
                
            end

        end % observations  
        
        gfp(r, :) =  std(squeeze(mean(data_rand, 1)), 1);
        
    end % randomnizations 
    
    % p-value for each feature
    pval = sum(gfp >= repmat(gfp(1, :), n_rand + 1, 1)) ./ n_rand;
     
    % Acceptable false discovery rate 
    % (acceptable fraction of false positives)
    q = 0.05;
    
    % Control the false discovery rate for 
    % a family of hypothesis 
    [decision, ~, ~, ~] = fdr_bh(pval, q, 'pdep', 'no');

end