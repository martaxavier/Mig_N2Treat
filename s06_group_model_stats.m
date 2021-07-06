% Import report APIs 
import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
        
% Create a title for the report
if flag.report ~= 0
    
    my_title = strcat(upper(reg_model), ' GROUP MODEL STATS');
    H1 = get_report_heading(1,my_title);
    add(R,H1)  
    
end

n_subjects = length(subjects);
n_metrics = length(metrics);

subj_stats = cell(n_metrics,1); 
                                
% ------------------------------------------------------------
% Retreive subjects models
% ------------------------------------------------------------

% Go through metrics 
for m = 1 : n_metrics

    metric = metrics(m);
    
    % Get parameters for 
    % current metric
    get_metric_pars; 
    
%     if strcmp(reg_model,'l21_1') && n_bands == 1
%         continue
%     end
    
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
        model_in = strcat(metric,'_','model', ...
            '_',cv_method,'.mat');
        load(fullfile(path_data_in(s,r),model_in));
        efp(:,s) = model.efp(2:end);    
        
    end % finish looping through subjects 
    
    subj_stats{m,1}= efp;
        
end % finish looping through metrics 

% ------------------------------------------------------------
% Compute group stats 
% ------------------------------------------------------------ 

for m = 1 : n_metrics
    
    metric = metrics(m);
    
%     if strcmp(reg_model,'l21_1') && n_bands == 1
%         continue
%     end
    
    % Broadcast the current pipeline stage 
    disp(strcat('Creating report of model group ', ...
        ' results for metric'," ",metric," ..."));
    
    % Get parameters for 
    % current metric
    get_metric_pars;
    
    % Specify directory where images are to be saved 
    path_img_metric_out = strcat(path_img_out(r),'\',metric);

    % Create directory where results are to be saved 
    if ~exist(path_img_metric_out, 'dir')
        mkdir(path_img_metric_out); 
    end    

    % Create directory where results are to be saved 
    if ~exist(path_data_out(r), 'dir')
        mkdir(path_data_out(r)); 
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
    save(fullfile(path_data_out(r),data_out),'gstats');  

    % Perform the topographic consistency test
    [prob_consistency, pval_consistency, ~] = ...
        test_topo_consistency...
        (cell2mat(subj_stats(m,1)), metric, n_rand, path_pars);
    
        
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
                    
        report_group_stats(group_stats,...
            thresh_model, metric, R, prob_consistency, ...
            pval_consistency, flag.report, path_img_metric_out, path_pars);
        
        report_group_models(cell2mat(subj_stats(m,1)), ...
            metric, R, prob_consistency, pval_consistency, ...
            path_img_metric_out, path_pars)

    end
    
end % finish looping through metrics 


%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% [group_stats] = perform_group_analysis(subj_efp)           
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