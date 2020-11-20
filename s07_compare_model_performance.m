
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
    
    % Perform ANOVA of the parameters specified in the 
    % variables 'model_fields' and 'optimal_fields'
    [box,mcomp] = perform_anova1(models, ...
        model_fields,path_data_out(r));
    [box_opt,mcomp_opt] = perform_anova1(optimals, ...
        optimal_fields,path_data_out(r));
    
    % Merge the two structs together 
    for f = 1 : length(optimal_fields)
        box.(optimal_fields(f)) = ...
            box_opt.(optimal_fields(f));
        mcomp.(optimal_fields(f)) = ...
            mcomp_opt.(optimal_fields(f));
    end
    
    if flag.report ~= 0
        
        % Create a title for this segment
        % of the report 
        my_title = strcat('1-WAY ANOVA -', ...
            " ",upper(reg_model));
        H1 = get_report_heading(1,my_title);
        add(R,H1)  

        % Plot and report ANOVA results 
        fields = [model_fields optimal_fields];
        report_anova(box,mcomp,fields,metrics, ...
            path_data_out,path_img_out);

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

% Perform 2-way ANOVA of the parameters specified in
% the variables 'model_fields' and 'optimal_fields'
[box2,mcomp2] = perform_anova2(models, ...
    model_fields,path_data_out(3));
[box2_opt,mcomp2_opt] = perform_anova2(optimals, ...
    optimal_fields,path_data_out(3));

% Merge the two structs together 
for f = 1 : length(optimal_fields)
    box2.(optimal_fields(f)) = ...
        box2_opt.(optimal_fields(f));
    mcomp2.(optimal_fields(f)) = ...
        mcomp2_opt.(optimal_fields(f));
end
    
if flag.report ~= 0

    % Create a title for this segment
    % of the report 
    my_title = '2-WAY ANOVA';
    H1 = get_report_heading(1,my_title);
    add(R,H1)  

    % Plot and report 2-way ANOVA results 
    report_anova(box2,mcomp2,fields,reg_models, ...
        path_data_out,path_img_out)
        
end
    
%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% perform_anova1(models,fields,path_data_out)             
% ============================================================

function [box,mcomp] = perform_anova1(models,fields,path_data_out)

%   [box,mcomp] = perform_anova1(models,fields,path_data_out) performs
%                 analysis of variance (ANOVA) on the 'model' fields 
%                 specified in 'fields'
%
%   INPUTS:
%
%       models          the struct containing the final model/the 
%                       model at each fold, for each subject and 
%                       for each of the conditions (groups) under 
%                   	comparison
%       fields      	the fields (parameters) in the struct 
%                   	'models' to which ANOVA is to be applied  
%       path_data_out 	the directory where results are to be saved 
%
%   OUTPUTS:
%
%       box             struct containing the values of each field, 
%                       for each cv fold and subject (rows) and for 
%                       each group (columns)
%       mcomp           struct containing, for each field, the pairwise  
%                       multiple comparison of the ANOVA results for  
%                       each pair of groups 
%   
    
% -------------------------------------------------
% Read input information  
% -------------------------------------------------

% Read size of the problem 
n_subjects = size(models,1);
n_groups = size(models,2);
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
    
    % Compute pairwise results of the multiple comparison
    % test - obtain p-values for the hypothesis test that
    % the corresponding pairwise mean difference is not 0
    % Cols 1 and 2 contain the indices of the two samples
    % being compared; col 3 is the lower confidence interval,
    % col 4 is the estimate, col 5 is the upper confidence 
    % interval, col 6 is the p-value 
    multcomp = multcompare(field_stats); 

    % Save multiple comparison results for current 
    % field in output path 
    data_out = strcat(field,'_multcomp.mat');
    save(fullfile(path_data_out,data_out),'multcomp');
    
    % Add the current field's values to the struct
    % 'box' + add the current field's pairwise 
    % comparison results to the struct 'mcomp'
    box.(field) = field_values;
    mcomp.(field) = multcomp;
    
end

end

% ============================================================
% report_anova(box,mcomp,fields,groups,path_data_out,path_img_out)             
% ============================================================

function [] = report_anova(box,mcomp,fields,groups, ...
    path_data_out,path_img_out)

% -------------------------------------------------
% Plot and save boxplots 
% -------------------------------------------------

n_fields = length(fields);
n_groups = length(groups);

for f = 1 : n_fields
    
    field = fields(f);
    multcomp = mcomp.(field); 
        
    % Create report heading for this section 
    my_title = upper(field);
    H2 = get_report_heading(2,my_title);
    add(R,H2);

    % Pre-allocate array of pairwise pvalues 
    pvalues = zeros(n_groups*n_groups,1);
    
    % Assign the row and column subscript that correspond 
    % to each linear index of the 'multcomp' matrix 
    rows = multcomp(:,1); cols = multcomp(:,2);
    
    % Convert the linear indices of 'multcomp' to subscript values
    % and assign the pvalues contained in the 6th column of 'multcomp'
    % to the new pvalues array 
    pvalues(sub2ind([n_groups n_groups],rows,cols)) = multcomp(:,6);
    
    % Reshape the pvalues array into an upper triangular matrix 
    pvalues = reshape(pvalues,[n_groups n_groups]);
         
    % Re-organize the pvalues into a table format
    table = array2table(pvalues,'VariableNames', ...
        upper(groups),'RowNames',upper(groups));
    
    % Save the pvalues table for the current field 
    % in the specified output path 
    data_out = strcat(field(f),'_pvalues.txt');
    writetable(table,fullfile(path_data_out,data_out), ...
        'WriteVariableNames',true,'WriteRowNames', ...
        true,'Delimiter','\t'); 
   
    % Create heading for the current section 
    my_title = strcat('P-values of the', ...
        ' Parameters Pairwise Comparison');
    H3 = get_report_heading(3,my_title);
    add(R,H3);
    
    % Add table of p-values to the report 
    T = Table(table);
    T.TableEntriesStyle = {FontFamily('Arial'), ...
        Width('1in'),Color('black'),Bold};
    T.TableEntriesVAlign = 'middle';T.HAlign = 'center';
    T.Style = {RowHeight('0.5in')}; add(R,T)
  
    % Create heading for the current section 
    my_title = strcat('Boxplot of the', ...
        ' Parameters Values');
    H3 = get_report_heading(3,my_title);
    add(R,H3);
    
    % Create a boxplot image for the current field 
    % And add image to the report 
    b = boxplot(box.(field),'Colors',[0 0.4470 0.7410]);
    set(b,{'linew'},{2}); 
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
    set(gca,'XTickLabel',labels_in,'FontSize',18)
    ylabel(upper(field),'FontSize',22);
    
    % Save the boxplot image in the specified output path
    img_out = strcat(field(f),'_boxplot.png');
    source = fullfile(path_img_out,img_out);
    saveas(gcf,source); I = Image(source);
    I.Style={ScaleToFit(true),HAlign('center')};
    add(R,I);
    
end

end

% ============================================================
% perform_anova1(models,fields,path_data_out)             
% ============================================================

function [box,mcomp] = perform_anova2(models,fields,path_data_out)

%   [box,mcomp] = perform_anova1(models,fields,path_data_out) performs
%                 analysis of variance (ANOVA) on the 'model' fields 
%                 specified in 'fields'
%
%   INPUTS:
%
%       models          the struct containing the final model/the 
%                       model at each fold, for each subject and 
%                       for each pair of conditions (groups) under 
%                   	comparison
%       fields      	the fields (parameters) in the struct 
%                   	'models' to which ANOVA is to be applied  
%       path_data_out 	the directory where results are to be saved 
%
%   OUTPUTS:
%
%       box             struct containing the values of each field, 
%                       for each cv fold and subject (rows) and for 
%                       each combination of groups (columns/depth)
%       mcomp           struct containing, for each field, the pairwise  
%                       multiple comparison of the ANOVA results for  
%                       each pair of groups 
%   
    
% -------------------------------------------------
% Read input information  
% -------------------------------------------------

% Read size of the problem 
n_subjects = size(models,1);
n_groups = [size(models,2) size(models,3)];
n_folds = length(models(1,1).(fields(1)));
n_fields = length(fields);

% Go through fields 
for f = 1 : n_fields
    
    field = fields(f);
    field_values = zeros(n_folds,...
        n_subjects,n_groups(1),n_groups(2));
    
    % Go through groups
    for g1 = 1 : n_groups(1)
        
        for g2 = 1 : n_groups(2)
        
            % Go through subjects 
            for s = 1 : n_subjects

                field_values(:,s,g1,g2) = ...
                    models(s,g1,g2).(field);

            end   
            
        end
        
    end
    
    % Prepare field values matrix for ANOVA 
    field_values = reshape(field_values, ...
        [n_folds*n_subjects*n_groups(1),n_groups(2)]);
    
    % Perform 2-way ANOVA for current field 
    reps = n_folds*n_subjects;
    [~,~,field_stats] = anova2(field_values,reps);
    
    % Compute pairwise results of the multiple comparison
    % test - obtain p-values for the hypothesis test that
    % the corresponding pairwise mean difference is not 0
    % Cols 1 and 2 contain the indices of the two samples
    % being compared; col 3 is the lower confidence interval,
    % col 4 is the estimate, col 5 is the upper confidence 
    % interval, col 6 is the p-value 
    multcomp = multcompare(field_stats); 

    % Save multiple comparison results for current 
    % field in output path 
    data_out = strcat(field,'_multcomp.mat');
    save(fullfile(path_data_out,data_out),'multcomp');
    
    % Add the current field's values to the struct
    % 'box' + add the current field's pairwise 
    % comparison results to the struct 'mcomp'
    box.(field) = field_values;
    mcomp.(field) = multcomp;
    
end

end
