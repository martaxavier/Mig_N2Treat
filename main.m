% Define and excecute analysis pipeline 

close all 
clear all 

%------------------------------------------------------------------
% Define pipeline steps 
%------------------------------------------------------------------

% LEVEL 0
flag.process_eeg = 0;
flag.process_bold = 0;
flag.process_bold_imgs = 0;
flag.extract_eeg_markers = 0;

% LEVEL 1 
flag.compute_features = 0;
flag.deconvolve_bold = 0;

% LEVEL 2
flag.frequency_analysis = 0;
flag.correlation_analysis = 0;

% LEVEL 3
flag.group_correlation_analysis = 0;

% LEVEL 4
flag.estimate_acf_order = 0;
flag.optimize_cv_pars = 0;

% LEVEL 5
flag.fit_models = 0;
flag.report_models = 0;

% LEVEL 6
flag.group_model_stats = 1;

% LEVEL 7
flag.compare_model_performance = 0;

% REPORT 
flag.report = 1;    % 2 to generate report + report images
                    % 1 to generate files + report + all images (slow)
                    % 0 to generate only output files 

%------------------------------------------------------------------
% Execute analysis pipeline 
%------------------------------------------------------------------

% Run configuration script
config

% Import report APIs 
if flag.report ~= 0
    
    % These statements eliminate the need to qualify
    % the names of DOM and REPORT objects, e.g., you 
    % can refer to mlreportgen.dom.Document as Document
    import mlreportgen.dom.*;
    import mlreportgen.report.*;
    iptgetpref('ImshowInitialMagnification');
   
    % Create report output directory 
    % if non existent 
    if ~exist(path.report, 'dir')
        mkdir(path.report)
    end
    
end

% Create parameters output 
% directory if non existent 
if ~exist(path.pars,'dir')
    mkdir(path.pars);
end

% -------------------------------------------------
% Process EEG
% -------------------------------------------------
if flag.process_eeg
    
    % I HAVE TO REMOVE THESE RETURNS
    % SINCE THEY LEAVE THE SCRIPT
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Define input/output data/path
    data_in = filename.eeg_raw; 
    data_out = filename.eeg_processed;
    path_data_in = path.eeg_raw;
    path_data_out = path.eeg_processed;
   
    % Go through subjects  
    s00_process_EEG;
    
end

% -------------------------------------------------
% Process BOLD
% -------------------------------------------------
if flag.process_bold
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Define input/output data/path
    data_in = filename.bold_raw; 
    data_out = filename.bold_processed;
    path_data_in = path.bold_raw;
    path_data_out = path.bold_processed;
    path_img_out = strcat('IMAGES\',path_data_out);
   
    % Run script 
    s00_process_BOLD;
    
end

% -------------------------------------------------
% Process BOLD images 
% -------------------------------------------------
if flag.process_bold_imgs
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Leave if non of the metrics
    % specified are performed with 
    % the deconvolved BOLD signal 
    if ~contains(metrics,'deconv')
        return
    end
    
    % Define input/output data/path
    data_bold_in = filename.bold_img_raw; 
    data_bold_out = filename.bold_img_processed;
    data_dmn_in = filename.bold_mask_dmn;
    data_dmn_out = filename.bold_mask_dmn_bin;
    path_data_in = path.bold_img_raw;
    path_data_out = path.bold_img_processed;
    path_img_out = strcat('IMAGES\',path_data_out);
   
    % Run script 
    s00_process_BOLD_imgs;
   
end

% -------------------------------------------------
% Extract EEG markers 
% -------------------------------------------------

if flag.extract_eeg_markers
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
               
    % Define input/output data/paths
    data_in = filename.eeg_processed;
    data_eeg_out = filename.eeg_markers;
    data_eeg_sub_task_out = filename.eeg_markers_sub_task;
    data_bold_out = filename.bold_markers;
    path_data_in = path.eeg_processed;
    path_data_eeg_out = path.eeg_markers;
    path_data_bold_out = path.bold_markers;
    
    % Run script
    s00_extract_EEG_markers;
        
end

% -------------------------------------------------
% Compute features 
% -------------------------------------------------
all_metrics = metrics;
if flag.compute_features
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Define input/output data/paths
    data_in = filename.eeg_processed;        
    markers_in = filename.eeg_markers;
    markers_sub_task_in = filename.eeg_markers_sub_task;
    path_data_in = path.eeg_processed;
    path_markers_in = path.eeg_markers;

    data_out = filename.eeg_feature;
    data_out_eeg_fs = filename.eeg_feature_eeg_fs;
    data_out_conv = filename.eeg_feature_conv;
    data_out_delay = filename.eeg_feature_delay;
    path_data_out = path.eeg_feature;
    path_img_out = strcat('IMAGES\',path_data_out,'\features');
        
    % Compute power features if any supported
    % power feature metric exists in metrics 
    if max(contains(power_metrics,all_metrics))
        
        % Find metrics that belong to the power_metrics 
        metrics = power_metrics(contains(power_metrics, ...
            all_metrics));

        % Run script
        s01_compute_power_features;
        
    end
    
    % Compute power features if any supported
    % connectivity feature metric exists in metrics 
    if max(contains(connectivity_metrics,all_metrics))
        
        % Find metrics that belong to the power_metrics 
        metrics = connectivity_metrics(contains(...
            connectivity_metrics,all_metrics));
        
        % Run script
        s01_compute_connectivity_features;
        
    end
    
end

metrics = all_metrics;

% -------------------------------------------------
% Deconvolve BOLD signal 
% -------------------------------------------------

if flag.deconvolve_bold
    
    % Leave if non of the metrics
    % specified are performed with 
    % the deconvolved BOLD signal 
    if ~contains(metrics,'deconv')
        return
    end
    
    % Assign input/output data/paths, according
    % to the BOLD deconvolution method defined 
    if strcmp(deconv_method,'time_series')
        
        data_in = filename.bold_processed;
        data_out = filename.bold_deconv;
        path_data_in = path.bold_processed;
        path_data_out = path.bold_deconv;
        
    elseif strcmp(deconv_method,'voxel_wise')
        
        data_in = filename.bold_img_processed;
        struct_in = filename.bold_img_raw;
        data_out = filename.bold_img_deconv;
        path_data_in = path.bold_img_processed;
        path_struct_in = path.bold_img_raw;
        path_data_out = path.bold_img_deconv;
        
    end
    
    dmn_in = filename.bold_mask_dmn_bin;
        
    path_img_out = strcat('IMAGES\',path_data_out);
    path_report_out = path.report;
        
    % Create the report object 
    if flag.report == 2
       
        my_report = 'RESULTS - BOLD DECONVOLUTION';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
    
    s01_deconvolve_bold;
    
    if flag.report == 2
       
        close(R) 
        
    end

end

% -------------------------------------------------
% Frequency analysis  
% -------------------------------------------------

if flag.frequency_analysis
    
    % Define input/output data/paths
    data_in = filename.eeg_feature_eeg_fs;
    path_data_in = path.eeg_feature;
    path_data_out = path.frequency;
    path_img_in = strcat('IMAGES\',path_data_in);
    path_img_out = strcat('IMAGES\',path_data_out);
    path_report_out = path.report;

    % Leave if image generation 
    % flag is turned off 
    if flag.report == 0
        
        return
     
    % Create report object 
    else
        
        my_report = 'RESULTS - FREQUENCY ANALYSIS';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
    
    % Run script
    s02_frequency_analysis;
    close(R);
            
end

% -------------------------------------------------
% Correlation analysis  
% -------------------------------------------------

if flag.correlation_analysis
    
    % Define input/output data/paths
    data_out = filename.correlation;
    path_eeg_in = path.eeg_feature;
    path_bold_in = path.bold_processed;
    path_data_out = path.correlation;
    path_img_out = strcat('IMAGES\',path_data_out);
    path_pars_in = path.pars;
    path_report_out = path.report;

    % Create report object
    if flag.report ~= 0
        
        my_report = 'RESULTS - CORRELATION ANALYSIS';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
    
    % Run script
    s02_correlation_analysis;
    
    if flag.report ~= 0
        close(R);
    end
            
end


if flag.group_correlation_analysis

    % Define input/output data/paths 
    data_in = filename.correlation;
    data_out = filename.correlation_group;
    path_data_in = path.correlation;
    path_data_out = path.correlation_group;
    path_img_out =  strcat('IMAGES\',path_data_out);
    path_report_out = path.report;
    
    % Create report object
    if flag.report ~= 0
        
        my_report = strcat('RESULTS -', ...
            ' GROUP CORRELATION ANALYSIS');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end

    % Create output directories if non existent 
    if ~exist(path_data_out, 'dir')
        mkdir(path_data_out)
    end
    
    if ~exist(path_img_out,'dir')
        mkdir(path_img_out)
    end

    s03_group_correlation_analysis;
    
    if flag.report ~=0
        close(R);
    end
    
end

% -------------------------------------------------
% Estimate order of the ACF model 
% -------------------------------------------------
if flag.estimate_acf_order
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Define input/output data/paths 
    data_in = filename.bold_processed;
    data_deconv_in = filename.bold_deconv;
    path_data_in = path.bold_processed;
    path_deconv_in = path.bold_deconv;
    path_pars_out = path.pars;
    path_img_out = strcat('\IMAGES\', ...
        path_pars_out);
    
    % Run script 
    s04_estimate_acf_order;
    
end

% -------------------------------------------------
% Optimize Cross-Validation Parameters 
% -------------------------------------------------
if flag.optimize_cv_pars
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Go through CV method 
    for r = 1 : length(reg_models)
        
        reg_model = reg_models(r);
        
        % Define input/output data/paths 
        path_eeg_in = path.eeg_feature;
        path_bold_in = path.bold_processed; 
        path_pars_in = path.pars;
        path_pars_out = path.pars;
                
        % Optimize CV parameters 
        % for current CV method 
        s04_optimize_cv_pars;
        
    end
    
end
    

% -------------------------------------------------
% Fit EEG-BOLD models
% -------------------------------------------------
if flag.fit_models
    
    % Leave if file generation 
    % flag is turned off 
    if flag.report == 2
        return
    end
    
    % Define input/output data/paths 
    path_eeg_in = path.eeg_feature;
    path_bold_in = path.bold_processed;
    path_data_out = path.model; 
    path_img_out = strcat('IMAGES\',path_data_out);
    path_pars_in = path.pars;
    path_report_out = path.report;

    % Go through CV methods
    for r = 1 : length(reg_models)
        
        reg_model = reg_models(r);

        % Fit model for current 
        % CV method 
        s05_fit_models;
        
    end
    
end

% -------------------------------------------------
% Report EEG-BOLD Models
% -------------------------------------------------
if flag.report_models
    
    % Leave if image generation
    % flag is turned off 
    if flag.report == 0
        return
    end 
    
    % Define input/output data/paths  
    path_eeg_in = path.eeg_feature;
    path_bold_in = path.bold_processed;
    path_data_in = path.model;
    
    path_img_out = strcat('IMAGES\',path_data_in);
    path_pars_in = path.pars;
    path_report_out = path.report;
    
    % Go through CV methods 
    for r = 1 : length(reg_models)
        
        reg_model = reg_models(r);

        % Create the report object
        my_report = strcat('RESULTS -', ...
            " ",upper(reg_model),' MODELS');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
       
        % Run script 
        s05_report_models;

        close(R)

    end 
    
end

% -------------------------------------------------
% Group Statistics of EEG-BOLD Models 
% -------------------------------------------------

if flag.group_model_stats
    
    % Define input/output data/paths  
    data_out = filename.model_group;
    path_data_in = path.model;
    path_data_out = path.model_group; 
    
    path_img_out = strcat('IMAGES\',path_data_out);
    path_pars_in = path.pars;
    path_report_out = path.report;
    
    % Go through CV methods 
    for r = 1 : length(reg_models)
        
        reg_model = reg_models(r);

        if flag.report ~= 0
            
            % Create the report object 
            my_report = strcat('RESULTS -', ...
                " ",upper(reg_model), ...
                ' GROUP MODELS');
            R = Report(my_report,'pdf');
            R.Layout.Landscape = true;
            R.OutputPath = strcat(path_report_out, ...
                '\',my_report);
            open(R)
  
        end
        
        % Run script 
        s06_group_model_stats;
        
        if flag.report ~= 0
            close(R);
        end
        
    end
    
end

% -------------------------------------------------
% Performance of EEG-BOLD Models 
% -------------------------------------------------

if flag.compare_model_performance
    
    % Define input/output data/paths
    path_data_in = path.model; 
    path_data_out = path.compare_performance; 
    path_img_out = strcat('IMAGES\',path_data_out);
    path_report_out = path.report;
    
     if flag.report ~= 0
            
        % Create the report object 
        my_report = strcat('RESULTS - ', ...
            ' MODELS PERFORMANCE COMPARISON');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
  
     end
        
    s07_compare_model_performance;
    
    if flag.report ~= 0
        close(R);
    end
    
end

% -------------------------------------------------
% Final report 
% -------------------------------------------------

% should contain: group stats of correlation and models  
% (images and topographical significance test) + 
% model performance + analysis pipeline + analysis parameters 
% (the last two should all be in the config script)
