% Performs group-level correlation analysis, from the results
% of the subject-level correlation analysis, in which the Pearson's
% correlation between the EEG features and the corresponding BOLD
% signal was estimated 
% Plots the results, saves them in the IMAGES directory and  
% generates an automatic report 

% Import report APIs 
import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
        
% Create a title for the report
if flag.report ~= 0
    
    my_title = 'GROUP CORRELATION ANALYSIS';
    H1 = get_report_heading(1,my_title);
    add(R,H1)  
    
end

n_subjects = length(subjects);
n_sessions = length(sessions); 
n_metrics = length(metrics);

% ------------------------------------------------------------
% Retreive subjects correlation stats
% ------------------------------------------------------------    

subj_stats = cell(n_metrics, 2); % first column rho
                                 % second column pval
% Retreive subject-level
% correlation stats 
for m = 1 : n_metrics
   
    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);
       
    % Pre-allocate correlation stats variables
    rho = zeros(n_features, n_subjects, n_sessions);
    pval_consistency = zeros(n_features, n_subjects, n_sessions);
    
    for se = 1 : n_sessions
        
        for s = 1 : n_subjects

            subject = subjects(s);

            % Load subject-level correlation results 
            corr_in = strcat(metric, '_', data_in);
            load(fullfile(path_data_in(s, se), corr_in), 'stats');
            rho(:, s, se) = stats.rho; 
            pval_consistency(:, s, se) = stats.pval; 

        end % finish looping through subjects
        
    end % finish looping through sessions 
      
    subj_stats{m, 1}= reshape(rho, [n_features n_subjects*n_sessions]);
    subj_stats{m, 2} = reshape(pval_consistency, ...
        [n_features n_subjects*n_sessions]);
    
end % finish looping through metrics 

% ------------------------------------------------------------
% Perform group level analyses  
% ------------------------------------------------------------  
 
% Go through metrics 
for m = 1 : n_metrics

    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);  
        
    % Broadcast the current pipeline stage 
    disp(strcat('Performing group correlation', ...
        ' analysis for metric'," ", metric, ' ...'));
        
    % Add metric heading 
    % to the report 
    if flag.report ~=0
        
        my_title = upper(metric);
        H2 = get_report_heading(2, my_title);
        add(R, H2);
        
    end   
       
    % ------------------------------------------------------
    % Group-level correlation analysis  
    % ------------------------------------------------------
    
    % Perform group-level stats 
    % using the subj-level stats 
    gstats = perform_group_analysis...
        (cell2mat(subj_stats(m, 1)), ...
        thresh_fdr);
    group_stats = [];
    group_stats(:,1) = gstats.tstat;
    group_stats(:,2) = gstats.decision; 
    
    % Save results of the group correlation
    % analysis for the current metric 
    corr_out = strcat(metric, '_', data_out);
    save(fullfile(path_data_out, corr_out), 'gstats');  

    % ------------------------------------------------------
    % Topographic consistency test  
    % ------------------------------------------------------
    
    % Perform the topographic consistency test
    [~, pval_consistency, ~] = ...
        topo_consistency_test...
        (cell2mat(subj_stats(m, 1)), ...
        n_rand, metric, path_pars);
    
    % ------------------------------------------------------
    % Plot group-level and topographic consistency results   
    % ------------------------------------------------------    
    
    % Generate and save plots of the 
    % group-level correlation analysis 
    % results 
    if flag.report ~= 0
        
         report_group_stats(group_stats, ...
             thresh_tstat, thresh_fdr, pval_consistency, ...
             metric, R, flag.report, path_img_out, path_pars);
         
    end
       
end % finish looping through metrics

%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% [group_stats] = perform_group_analysis(subject_stats)           
% ============================================================

function [group_stats] = perform_group_analysis(subj_stats_rho, thresh_fdr)
%
%   [group_stats] = perform_group_stats performs second-level  
%                   analysis from the results of the first-level
%                   (subject-level) correlation statistics, using
%                   a one-sample t-test to test for significance 
%
%   Inputs:
%
%     subject_stats     first-level correlation statistics 
% 
%    Outputs:
%
%    group_stats        second-level correlation statistics 
%

[~, pval, ~, group_stats] = ttest(subj_stats_rho');

% Correct for multiple comparisons, using the FDR correction 
[decision, ~, ~, ~] = fdr_bh(pval, thresh_fdr, 'pdep', 'no');

group_stats.pval = pval;
group_stats.decision = decision; 

end