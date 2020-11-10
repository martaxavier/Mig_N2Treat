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
n_metrics = length(metrics);

% ------------------------------------------------------------
% Retreive subjects correlation stats
% ------------------------------------------------------------    

subj_stats = cell(n_metrics,2); % first column rho
                                % second column pval
% Retreive subject-level
% correlation stats 
for m = 1 : n_metrics
   
    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);
       
    % Pre-allocate correlation stats variables
    rho = zeros(n_features,n_subjects);
    pval = zeros(n_features,n_subjects);
    
    for s = 1 : n_subjects
        
        subject = subjects(s);
        
        % Load subject-level correlation results 
        corr_in = strcat(metric,'_',data_in);
        load(fullfile(path_data_in(s),corr_in),'stats');
        rho(:,s) = stats.rho; pval(:,s) = stats.pval; 
        
    end % finish looping through subjects
      
    subj_stats{m,1}= rho;
    subj_stats{m,2} = pval;
    
end % finish looping through metrics 

%-----------------------------------------------------------
% Go through metrics 
%-----------------------------------------------------------
                            
for m = 1 : n_metrics

    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);  
        
    % Broadcast the current pipeline stage 
    disp(strcat('Performing group correlation', ...
        ' analysis for metric'," ", metric, ' ...'));
        
    % Perform group-level stats 
    % using the subj-level stats 
    gstats = perform_group_analysis...
        (cell2mat(subj_stats(m,1)));
    group_stats = [];
    group_stats(:,1) = gstats.tstat;
    group_stats(:,2) = gstats.decision; 
    
    % Save results of the group correlation
    % analysis for the current metric 
    corr_out = strcat(metric,'_',data_out);
    save(fullfile(path_data_out,corr_out),'gstats');  

    % Perform the topographic consistency test
    [prob,but_one_prob] = ...
        test_topo_consistency...
        (cell2mat(subj_stats(m,1)),metric,n_rand);

    % Add metric heading 
    % to the report 
    if flag.report ~=0
        
        my_title = upper(metric);
        H2 = get_report_heading(2,my_title);
        add(R,H2);
        
    end

    % Generate and save plots of the 
    % group-level correlation analysis 
    % results 
    if flag.report ~= 0
        
        plot_group_stats(group_stats, ...
            thresh_corr,metric,R,prob,flag.report, ...
            path_img_out);
        
    end

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

[decision,pval,~,group_stats] = ttest(subj_stats_rho');
group_stats.pval = pval;
group_stats.decision = decision; 

end