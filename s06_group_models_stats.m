% Import report APIs 
import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
        
% Create a title for the report
if flag.report ~= 0
    
    my_title = strcat(model, ' GROUP MODEL STATS');
    H1 = get_report_heading(1,my_title);
    add(R,H1)  
    
end

n_subjects = length(subjects);
n_metrics = length(metrics);

subj_stats = cell(n_metrics,1); 
                                
% ------------------------------------------------------------
% Retreive subjects models
% ------------------------------------------------------------

for m = 1 : n_metrics

    metric = metrics(m);
    
    % Get parameters for 
    % current metric
    get_metric_pars; 
    
    % Compute # of features
    n_features = prod(dim);
    
    % Pre-allocate EFP variable 
    efp = zeros(n_features,n_subjects);
    
    %--------------------------------------------------------    
    % Go through subjects 
    %--------------------------------------------------------
                   
    for s = 1 : n_subjects
    
        subject = subjects(s);
        
        % Load model results for current subject,
        % model pair 
        model_in = strcat(metric,'_','model.mat');
        load(fullfile(path_model_in(s,m),model_in));
        efp(:,s) = model.efp(2:end);    
        
    end % finish looping through subjects 
    
    subj_stats{m,1}= model;
        
end % finish looping through metrics 

% ------------------------------------------------------------
% Compute group stats 
% ------------------------------------------------------------ 

for m = 1 : n_metrics
    
    metric = metrics(m);
    
    % Broadcast the current pipeline stage 
    disp(strcat('Creating report of model group ', ...
        ' results for metric'," ",metric," ..."));
    
    % Get parameters for 
    % current metric
    get_metric_pars;
    
    % Specify directory where images are to be saved 
    path_img_metric_out = strcat(path_img_out,'\',metric);

    % Create directory where results are to be saved 
    if ~exist(path_img_metric_out, 'dir')
        mkdir(path_img_metric_out); 
    end    

    % Perform group-level stats using
    % the subject-level stats 
    gstats = perform_group_analysis...
        (cell2mat(subj_stats(m,1)));
    group_stats = [];
    group_stats(:,1) = gstats.tstat;
    group_stats(:,2) = gstats.decision; 
    
    % Save results of the group correlation
    % analysis for the current metric 
    model_out = strcat(metric,'_',data_out);
    save(fullfile(path_data_out,data_out),'gstats');  

    % Perform the topographic consistency test
    [prob,but_one_prob] = ...
        test_topo_consistency...
        (cell2mat(subj_stats(m,1)),metric,n_rand);
    
        
    % Add metric heading
    % to the report 
    if flag.report ~= 0
        
        my_title = upper(metric);
        H2 = get_report_heading(2,my_title);
        add(R,H2);
        
    end
    
    % Generate and save plots of the 
    % group-level model statistics 
    % results 
    if flag.report ~= 0
        
        plot_group_stats(group_stats,...
            thresh_model,metric,R,prob, ...
            flag.report,path_img_out);
        
    end
    
end % finish looping through metrics 


%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% [group_stats] = perform_group_analysis(subject_stats)           
% ============================================================

function [group_stats] = perform_group_analysis(subj_efp)
%
%   [group_stats] = perform_group_stats performs second-level  
%                   analysis from the results of the first-level
%                   (subject-level) model statistics, using
%                   a one-sample t-test to test for significance 
%
%   Inputs:
%
%     subj_efp          first-level regression statistics 
% 
%    Outputs:
%
%    group_stats        second-level regression statistics 
%

[decision,pval,~,group_stats] = ttest(subj_efp');
group_stats.pval = pval;
group_stats.decision = decision; 

group_stats.pval(isnan(group_stats.tstat))=0;
group_stats.decision(isnan(group_stats.tstat))=0;
group_stats.tstat(isnan(group_stats.tstat))=0;

end