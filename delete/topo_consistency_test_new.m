function [prob, pval, but_one_prob] = topo_consistency_test_new(...
    first_stats, metric, n_rand, path_pars)
%
%   [prob] = topo_consistency_test tests, for the current metric, the 
%            consistency of the group-level topographies, through  
%            a global randomization test 
%                 
%
%   Inputs:
%
%     first_stats	    the first-level correlation/model values
%     metric                the current EEG metric 
%     n_rand                the number of randomizations 
% 
%    Outputs:
%
%    prob                   the probability of rejecting the 
%                           null-hypothesis 
%    prob_but_one           the probability of rejecting the 
%                           null-hypothesis, after the removal
%                           of each element of the group 
%

% Get metric parameters 
% for the current metric
get_metric_pars;

% First-level stats 
rho = first_stats';

% Number of elements in the group
n_group = size(rho, 1);

% Reshape first-level stats according 
% to the current metric's specified 
% feature space dimensions 
rho = reshape(rho, [n_group dim]);

% Generate #n_rand new rho matrices,
% the same size as the input matrix,
% for each element of the group 
[prob, pval] = global_rand_test(rho, n_rand, dim);

but_one_prob = zeros([n_subjects size(prob)]);

% Test the impact that removing each
% of the group elements  has on the 
% estimated probability of significance 
% for g = 1 : n_group
%    group_elem = 1 : n_group;
%    group_elem(g) = [];
%    but_one_rho = ...
%        squeeze(rho(group_elem,:,:,:,:));
%    but_one_prob(s,:,:) = ...
%        global_rand_test(but_one_rho, n_rand, dim);
% end

end

% ============================================================
% [rand_rho] = global_rand_test(rho, n_rand, dim)         
% ============================================================

function [prob, pval] = global_rand_test(rho, n_rand, dim)

    n_group = size(rho,1);
    n_chans = dim(1);
    n_delays = dim(2);
    n_bands = dim(3);

    rand_rho = zeros([n_group dim n_rand]);

    for r = 1 : n_rand  
        for g = 1 : n_group
            for ii = 1 : n_delays*n_bands
                p = randperm(n_chans);
                [band, delay] = ...
                    ind2sub([n_delays n_bands], ii);
                rand_rho(g, :, delay, band, r) = rho(g, p, band, delay);
            end % topographies 
        end % group elements
    end % randomnizations
    
    all_rho = cat(5, rho, rand_rho);

    % For each delay-band pair, compute
    % the average correlation values of 
    % the group
    mean_rho = squeeze(mean(all_rho));

    % Compute the root mean square of 
    % the average correlation values
    rms_rho = squeeze(rms(mean_rho));
    
    % Compute the global field power of 
    % the average correlation values 
    %gfp_rho = squeeze(std(mean_rho));

    prob = 100*sum( rms_rho(:, :, 2:end)...
        - rms_rho(:, :, 1) > 0, 3) / n_rand;
    
    is_greater = rand_rho > repmat(rho, ...
        1, 1, 1, 1, n_rand);
    pval = sum(is_greater, 5) ./ n_rand;

end