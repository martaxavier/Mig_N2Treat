
% ============================================================
% [prob,but_one_prob] = test_topo_consistency(subj_stats_rho, ...
%                       metric,n_rand)        
% ============================================================

function [prob,but_one_prob] = test_topo_consistency(...
    subj_stats,metric,n_rand)
%
%   [prob] = test_topo_consistency tests, for each metric, the 
%            consistency of the group-level topographies, through  
%            a global randomization test 
%                 
%
%   Inputs:
%
%     subject_stats         the subject-level correlation/model values
%     metric                the current EEG-BOLD metric
%     n_rand                the number of randomizations 
% 
%    Outputs:
%
%    prob                   the probability of rejecting the 
%                           null-hypothesis 
%    prob_but_one           the probability of rejecting the 
%                           null-hypothesis, after the removal
%                           of one subject (one value for each)
%

% Get metric parameters 
% for the current metric
get_metric_pars;

% Organize the input data 
% Rho represents the subject-level
% correlation or model values 
rho = subj_stats';
n_subjects = size(rho,1);
rho = reshape(rho,[n_subjects dim]);

% Generate #n_rand new rho matrices,
% the same size as the input matrix,
% for each subject 
prob = global_rand_test(rho,n_rand,dim);

but_one_prob = zeros([n_subjects size(prob)]);

% Test the impact that removing each
% of the subjects has on the estimated
% probability of significance 
% for s = 1 : n_subjects
%    subjs = 1:n_subjects;
%    subjs(s)=[];
%    but_one_rho = ...
%        squeeze(rho(subjs,:,:,:,:));
%    but_one_prob(s,:,:) = ...
%        global_rand_test(but_one_rho,n_rand,dim);
% end

end

% ============================================================
% [rand_rho] = global_rand_test(rho,n_rand,dim)         
% ============================================================

function [prob]=global_rand_test(rho,n_rand,dim)

    n_subjects = size(rho,1);
    n_chans = dim(1);
    n_delays = dim(2);
    n_bands = dim(3);

    rand_rho = zeros([n_subjects dim n_rand]);

    for r = 1 : n_rand  
        for s = 1 : n_subjects
            for ii = 1 : n_delays*n_bands
                p = randperm(n_chans);
                [row,col] = ...
                    ind2sub([n_delays n_bands],ii);
                rand_rho(s, :, row, col, r) = rho(s,p,row,col);
            end   
        end
    end
    
    all_rho = cat(5,rho,rand_rho);

    % For each delay-band pair, compute
    % the average correlation values of 
    % the group
    mean_rho = squeeze(mean(all_rho));

    % Compute the root mean square of 
    % the average correlation values
    rms_rho = squeeze(rms(mean_rho));

    prob = 100*sum(rms_rho(:,:,2:end)...
        -rms_rho(:,:,1)>0,3)/n_rand;

end