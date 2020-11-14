
% Compares the model performance between all of the metrics
% under analysis, for each of the regression models (1-way ANOVA)
% Compares the model performance between all of the metrics  
% and regression models under analysis (2-way ANOVA)

% Import report APIs 
import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 

% Specify parameters to be compared
% between conditions by ANOVA 
model_fields = "df";
optimal_fields = ["bic","nmse","corr"];

n_subjects = length(subjects);
n_metrics = length(metrics);
n_reg_models = length(reg_models);

% Load a random model and optimal
% object to retrieve fields 
model_in = strcat(metrics(1),'_', ...
    'model','_',cv_method,'.mat');
load(fulltile(path_model_in(1,1), ...
    model_in),'model','optimal');
            
% Pre-allocate arrays of models and optimals             
all_models = repmat(model, ...
    n_subjects,n_metrics,n_reg_models);
all_optimals = repmat(optimal, ...
    n_subjects,n_metrics,n_reg_models);

% n_opt_fields = length(opt_fields);
% opt = [opt_fields;ones(1,n_opt_fields)]; 
% opt = cellstr(opt(:)');
% all_optimals = struct(opt{:});
    
% Go through regression models
for r = 1 : n_reg_models
    
    % Go through metrics 
    for m = 1 : n_metrics

        metric = metrics(m);

        % Get parameters for 
        % current metric
        get_metric_pars; 

        %--------------------------------------------------------    
        % Go through subjects 
        %--------------------------------------------------------

        for s = 1 : n_subjects

            subject = subjects(s);

            % Load model results for current subject,
            % model pair 
            model_in = strcat(metric,'_','model', ...
                '_',cv_method,'.mat');
            model_folds_in = strcat(metric,'_', ...
                'model_folds_',cv_method,'.mat');             
            load(fullfile(path_model_in(s,r), ...
                model_in),'model');
            load(fullfile(path_model_in(s,r), ...
                model_folds_in),'optimal');
            
%             for f = 1 : n_opt_fields
%                 all_optimals.(opt_fields(f))(:,s,m,r) ...
%                     = optimal.(opt_fields(f));
%             end

            all_optimals(s,m,r) = optimal;
            all_models(s,m,r) = model;

        end % finish looping through subjects 

    end % finish looping through metrics 

end % finish looping through regression models

% ------------------------------------------------------------
% Perform 1-way ANOVA and report results 
% ------------------------------------------------------------

% Compare between different metrics 
% Go through regression models
for r = 1 : n_reg_models
    
    % Leave if there are not 
    % multiple metrics to compare 
    if n_metrics == 1
        return
    end
    
    reg_model = reg_models(r);
    models = squeeze(all_models(:,:,r));
    optimals = squeeze(all_optimals(:,:,r));
    
    % Perform ANOVA of the specified 
    % parameters in models 
    [box,mcmp] = perform_anova1(models, ...
        mod_fields,path_data_out);
    [box_cat,mcmp_cat] = perform_anova1...
        (optimals,optimal_fields,metrics,path_data_out);
    
    
    if flag.report ~= 0
        
        % Create a title for this segment
        % of the report 
        my_title = strcat('1-WAY ANOVA -', ...
            " ",reg_model);
        H1 = get_report_heading(1,my_title);
        add(R,H1)  

        % Plot and report ANOVA results 
        report_anova(box,mcomp,path_img_out);

    end
    
end

% ------------------------------------------------------------
% Perform 2-way ANOVA and report results 
% ------------------------------------------------------------

% Leave if there are not 
% multiple regression models 
if n_reg_models == 1
    return
end

varargin = {'metrics',metrics,'reg_models',reg_models};
anova_2 = perform_anova2(models,optimals,varargin, ...
    path_data_out);

if flag.report ~= 0

    % Create a title for this segment
    % of the report 
    my_title = '2-WAY ANOVA';
    H1 = get_report_heading(1,my_title);
    add(R,H1)  

    % Plot and report ANOVA results 
    report_anova(anova_2,path_img_out);

end
    
%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% perform_anova()             
% ============================================================

function [box,mcomp] = perform_anova1(models,fields,groups,path_data_out)

%   [] = perform_anova performs anova on the input data  
%        and plots and saves the corresponding boxplots 
%
%   INPUTS:
%

%     models           the struct containing the final model
%                      results for each subject and for each of
%                      the conditions being compared
%     optimals         the struct containing the model results 
%                      at each cv fold for each subject and for
%                      each of the conditions being compared
%     name             the name of the group being compared 
%     labels           labels of the conditions being compared 
%     path_imt_out     the directory where results are to be saved 
%   

% -------------------------------------------------
% Read input information  
% -------------------------------------------------

% Read size of the problem 
n_subjects = size(optimals,1);
n_groups = size(optimals,2);
n_folds = length(models(1,1).(fields(1)));
n_fields = length(fields);

% Go through fields 
for f = 1 : n_fields
    
    field = fields(f);
    field_values = zeros(n_folds,...
        n_subjects,n_groups);
    
    % Go through groups
    for g = 1 : n_groups
        
        % Go through subjects 
        for s = 1 : n_subjects
            
            field_values(:,s,g) = ...
                models(s,g).(field);
            
        end   
        
    end
    
    % Prepare field values matrix for ANOVA 
    field_values = reshape(field_values, ...
        [n_folds*n_subjects,n_groups]);
    
    % Perform ANOVA for current field 
    [~,~,field_stats] = anova1(field_values);
    
    % Compute pairwise results of the multiple
    % comparison test - obtain p-values for the
    % hypothesis test that the corresponding
    % pairwise mean difference is not equal to 0
    multcomp = multcompare(field_stats);
    
    % Save multiple comparison matrix for current 
    % field in output path 
    % WE SHOULD TURN THIS INTO A TABLE, WITH THE 
    % GROUP NAMES ABOVE - THIS COULD ALSO BE PRINTED
    % INTO THE REPORT 
    data_out = strcat(field(f),'_multcomp.mat');
    save(fullfile(path_data_out,data_out),'multcomp');
    
    % Add field values to box matrix 
    % for later plotting 
    box.(field(f)) = field_values;
    mcomp.(field(f)) = multcomp;
    
end

end

function [] = report_anova(fields,field_values,path_img_out)

% -------------------------------------------------
% Plot and save boxplots 
% -------------------------------------------------

% bic
b = boxplot(bic,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2}); 
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('BIC','FontSize',22);
img_out = strcat('anova_bic');
saveas(gcf,fullfile(path_out,img_out),'png');

% nmse
b = boxplot(nmse,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2});
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('NMSE','FontSize',22);
img_out = strcat('anova_nmse');
saveas(gcf,fullfile(path_out,img_out),'png');

% corr
b = boxplot(corr,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2});
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('CORR','FontSize',22);
img_out = strcat('anova_corr');
saveas(gcf,fullfile(path_out,img_out),'png');

% df
b = boxplot(df,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2});
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('DF','FontSize',22);
img_out = strcat('anova_df');
saveas(gcf,fullfile(path_out,img_out),'png');

end
