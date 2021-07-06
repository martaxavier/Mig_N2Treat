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
    
    my_title = 'GROUP CORRELATION ANALYSIS (SESSION)';
    H1 = get_report_heading(1, my_title);
    add(R, H1)  
    
end

n_sessions = length(sessions);
n_subjects = length(subjects);
n_metrics = length(metrics);

% ------------------------------------------------------------
% Retreive subjects correlation stats
% ------------------------------------------------------------    

se_stats = cell(n_metrics, n_subjects, 2); % first column rho
                                            % second column pval
% Retreive subject-level
% correlation stats 
for m = 1 : n_metrics
   
    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);
       
    for s = 1 : n_subjects
              
        % Create output directories if non existent 
        if ~exist(path_data_out(s), 'dir'); mkdir(path_data_out(s)); end
        if ~exist(path_img_out(s),'dir'); mkdir(path_img_out(s)); end
        
        % Pre-allocate correlation stats variables
        rho = zeros(n_features, n_sessions);
        pval_consistency = zeros(n_features, n_sessions);

        for se = 1 : n_sessions

            % Load subject-level correlation results 
            corr_in = strcat(metric, '_', data_in);
            load(fullfile(path_data_in(s, se), corr_in), 'stats');
            rho(:, se) = stats.rho; pval_consistency(:, se) = stats.pval; 

        end % finish looping through sessions
        
        se_stats{m, s, 1}= rho;
        se_stats{m, s, 2} = pval_consistency;        

    end % finish looping through subjects 
    
end % finish looping through metrics 

%-----------------------------------------------------------
% Go through metrics 
%-----------------------------------------------------------
                            
for m = 1 : n_metrics

    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);  
    
    % Add metric heading 
    % to the report 
    if flag.report ~=0

        my_title = upper(metric);
        H2 = get_report_heading(2, my_title);
        add(R, H2);

    end    
        
    for s = 1 : n_subjects
        
        subject = subjects(s);
        
        % Broadcast the current pipeline stage 
        disp(strcat('Performing session correlation', ...
            ' analysis for metric', " ", metric, ' subject', ...
            " ", subject, ' ...'));
        
        % Add subject heading 
        % to the report 
        if flag.report ~=0

            my_title = upper(subject);
            H3 = get_report_heading(3, my_title);
            add(R, H3);

        end           

        % Perform group-level stats 
        % using the subj-level stats 
        gstats = perform_group_analysis...
            (cell2mat(se_stats(m, s, 1)), ...
        thresh_fdr);
        group_stats = [];
        group_stats(:,1) = gstats.tstat;
        group_stats(:,2) = gstats.decision; 

        % Save results of the group correlation
        % analysis for the current metric 
        corr_out = strcat(metric, '_', data_out);
        save(fullfile(path_data_out(s), corr_out), 'gstats');  

        % Perform the topographic consistency test
        [~, pval_consistency, ~] = ...
            topo_consistency_test...
            (cell2mat(se_stats(m, s, 1)), n_rand, metric, path_pars);

        % Generate and save plots of the 
        % group-level correlation analysis 
        % results 
        if flag.report ~= 0

            report_group_stats(group_stats, ...
                thresh_tstat, thresh_fdr, pval_consistency, ...
                metric, R, flag.report, path_img_out(s), path_pars);

        end
        
    end % finish looping through subjects 

end % finish looping through metrics

%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% [group_stats] = perform_group_analysis(subject_stats)           
% ============================================================

function [group_stats] = perform_group_analysis(subj_stats_rho)
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