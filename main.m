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
flag.remove_ecg_save_chanlocs = 0;
flag.extract_eeg_markers = 0;

% LEVEL 1 
flag.compute_features = 0;
flag.deconvolve_bold = 0;

% LEVEL 2
flag.power_analysis = 0;
flag.connectivity_analysis = 0;
flag.correlation_analysis = 0;

% LEVEL 3
flag.group_correlation_analysis = 0;
flag.group_correlation_analysis_ses = 0; 
flag.power_reliability_analysis = 0;
flag.connectivity_reliability_analysis = 0;
flag.correlation_reliability_analysis = 1;

% LEVEL 4
flag.estimate_acf_order = 0;
flag.optimize_cv_pars = 0;

% LEVEL 5
flag.fit_models = 0;
flag.report_models = 0;

% LEVEL 6
flag.group_model_stats = 0;

% LEVEL 7
flag.compare_model_performance = 0;

% REPORT 
flag.report = 0;    % 2 to generate files + report + report images (fast,default)
                    % 1 to generate files + report + all images (slow)
                    % 0 to generate only output files (no report or images)

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
path_pars = path.pars;

% -------------------------------------------------
% Process EEG
% -------------------------------------------------
if flag.process_eeg
    
    % Define input/output data/path
    data_in = filename.eeg_raw; 
    data_out = filename.eeg_preproc;
    path_data_in = path.eeg_raw;
    path_data_out = path.eeg_preproc;
   
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
        
        s00_process_EEG;
   
    end
    
    fclose('all');
    
end

% -------------------------------------------------
% Process BOLD
% -------------------------------------------------
if flag.process_bold
    
    % Define input/output data/path
    data_in = filename.bold_raw; 
    data_out = filename.bold_preproc;
    path_data_in = path.bold_raw;
    path_data_out = path.bold_preproc;
    path_img_out = strcat('IMAGES\', path_data_out);
   
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
        
        % Run script 
        s00_process_BOLD;
    
    end % sessions 
    
end

% -------------------------------------------------
% Process BOLD images 
% -------------------------------------------------
if flag.process_bold_imgs
    
    % Define input/output data/path
    data_bold_in = filename.bold_img_raw; 
    data_bold_out = filename.bold_img_preproc;
    data_dmn_in = filename.bold_mask_dmn;
    data_dmn_out = filename.bold_mask_dmn_bin;
    path_data_in = path.bold_img_raw;
    path_data_out = path.bold_img_preproc;
    path_img_out = strcat('IMAGES\', path_data_out);
   
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
        
        % Run script 
        s00_process_BOLD_imgs;
        
    end % sessions 
   
end

% -------------------------------------------------
% Remove ECG and save chanlocs
% -------------------------------------------------

if flag.remove_ecg_save_chanlocs
               
    % Define input/output data/paths
    data_in = filename.eeg_preproc;
    data_out = data_in;
    path_data_in = path.eeg_preproc;
    path_data_out = path_data_in;
    
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
        
        % Run script
        s00_remove_ecg_save_chanlocs;
    
    end
        
end

% -------------------------------------------------
% Extract EEG markers 
% -------------------------------------------------

if flag.extract_eeg_markers
               
    % Define input/output data/paths
    data_in = filename.eeg_preproc;
    data_eeg_out = filename.eeg_markers;
    data_eeg_sub_task_out = filename.eeg_markers_sub_task;
    data_bold_out = filename.bold_markers;
    path_data_in = path.eeg_preproc;
    path_data_eeg_out = path.eeg_markers;
    path_data_bold_out = path.bold_markers;
    
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
    
        % Run script
        s00_extract_EEG_markers;
    
    end % sessions 
        
end

% -------------------------------------------------
% Compute features 
% -------------------------------------------------
all_metrics = metrics;
if flag.compute_features
    
    % Define input/output data/paths
    data_in = filename.eeg_preproc;        
    markers_in = filename.eeg_markers;
    markers_sub_task_in = filename.eeg_markers_sub_task;
    path_data_in = path.eeg_preproc;
    path_markers_in = path.eeg_markers;
    path_pars = path.pars;

    data_out = filename.eeg_feature;
    data_out_eeg_fs = filename.eeg_feature_eeg_fs;
    data_out_conv = filename.eeg_feature_conv;
    data_out_delay = filename.eeg_feature_delay;
    path_data_out = path.eeg_feature;
    path_img_out = strcat('IMAGES\', path_data_out);
        
    % Compute power features if any supported
    % power feature metric exists in metrics 
    if max(contains(power_metrics, all_metrics))
        
        % Find metrics that belong to the power_metrics 
        metrics = power_metrics(contains(power_metrics, ...
            all_metrics));

        % Go through sessions 
        for se = 1 : max(1, length(sessions))
        
            % Run script
            s01_compute_power_features;
        
        end % sessions 
        
    end
    
    % Compute power features if any supported
    % connectivity feature metric exists in metrics 
    if max(contains(connectivity_metrics,all_metrics))
        
        % Find metrics that belong to the power_metrics 
        metrics = connectivity_metrics(contains(...
            connectivity_metrics,all_metrics));
        
        % Go through sessions 
        for se = 1 : max(1, length(sessions))
        
            % Run script
            s01_compute_connectivity_features;
            
        end % sessions 
        
    end
    
end

metrics = all_metrics;

% -------------------------------------------------
% Deconvolve BOLD signal 
% -------------------------------------------------

if flag.deconvolve_bold
    
    % Assign input/output data/paths, according
    % to the BOLD deconvolution method defined 
    if strcmp(deconv_method,'time_series')
        
        data_in = filename.bold_preproc;
        data_out = filename.bold_deconv;
        path_data_in = path.bold_preproc;
        path_data_out = path.bold_deconv;
        
    elseif strcmp(deconv_method,'voxel_wise')
        
        data_in = filename.bold_img_preproc;
        struct_in = filename.bold_img_raw;
        data_out = filename.bold_img_deconv;
        path_data_in = path.bold_img_preproc;
        path_struct_in = path.bold_img_raw;
        path_data_out = path.bold_img_deconv;
        
    end
    
    dmn_in = filename.bold_mask_dmn_bin;
        
    path_img_out = strcat('IMAGES\', path_data_out);
    path_report_out = path.report;
        
    % Create the report object 
    % and open report 
    if flag.report == 2
       
        my_report = 'RESULTS - BOLD DECONVOLUTION';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
      
    % Run script 
    s01_deconvolve_bold;
    
    % Close report 
    if flag.report == 2       
        close(R) 
    end

end

% -------------------------------------------------
% Power analysis  
% -------------------------------------------------

if flag.power_analysis
    
    % NOTE: power analysis is performed 
    % in the frequency bands specified by  
    % the metric 
    
    % Define input/output data/paths
    data_in = filename.eeg_preproc;
    path_data_in = path.eeg_preproc;
    markers_in = filename.eeg_markers;
    path_markers_in = path.eeg_markers;      
    data_out = filename.power; 
    path_data_out = path.power;
    path_img_in = strcat('IMAGES\', path.eeg_feature);
    path_img_out = strcat('IMAGES\',path_data_out);
    path_report_out = path.report;

    % Create report object 
    % and open report 
    if flag.report ~= 0
        
        my_report = 'RESULTS - POWER ANALYSIS';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
    
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
        
        session = sessions(se);

        % Run script
        s02_power_analysis;
        
    end % sessions 
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
            
end

% -------------------------------------------------
% Connectivity analysis (static) 
% -------------------------------------------------

if flag.connectivity_analysis
    
    % NOTE: connectivity analysis is performed 
    % in the frequency bands specified by  
    % the metric 
    
    % Define input/output data/paths
    data_in = filename.eeg_preproc;
    path_data_in = path.eeg_preproc;
    markers_in = filename.eeg_markers;
    path_markers_in = path.eeg_markers;      
    data_out = filename.connectivity; 
    path_data_out = path.connectivity;
    path_img_in = strcat('IMAGES\', path.eeg_feature);
    path_img_out = strcat('IMAGES\',path_data_out);
    path_report_out = path.report;

    % Create report object 
    % and open report 
    if flag.report ~= 0
        
        my_report = 'RESULTS - CONNECTIVITY ANALYSIS';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
    
    % Go through sessions 
    for se = 1 : max(1, length(sessions))
        
        session = sessions(se);

        % Run script
        s02_connectivity_analysis;
        
    end % sessions 
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
            
end

% -------------------------------------------------
% Power reliability analysis 
% -------------------------------------------------

if flag.power_reliability_analysis

    % Define input/output data/paths 
    data_in = filename.power;
    data_out = filename.power_reliability;  
    path_data_in = path.power;
    path_data_out = path.power_reliability;
    path_img_out =  strcat('IMAGES\', path_data_out);
    path_report_out = path.report;
    
    % Create report object
    % and open report 
    if flag.report ~= 0
        
        my_report = strcat('RESULTS -', ...
            ' POWER RELIABILITY ANALYSIS');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end

    % Create output directories if non existent 
    if ~exist(path_data_out, 'dir'); mkdir(path_data_out); end
    if ~exist(path_img_out,'dir'); mkdir(path_img_out); end

    % Run script 
    s03_power_reliability_analysis;
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
    
end

% -------------------------------------------------
% Connectivity reliability analysis 
% -------------------------------------------------

if flag.connectivity_reliability_analysis

    % Define input/output data/paths 
    data_in = filename.connectivity;
    data_out = filename.connectivity_reliability;  
    path_data_in = path.connectivity;
    path_data_out = path.connectivity_reliability;
    path_img_out =  strcat('IMAGES\', path_data_out);
    path_report_out = path.report;
    
    % Create report object
    % and open report 
    if flag.report ~= 0
        
        my_report = strcat('RESULTS -', ...
            ' CONNECTIVITY RELIABILITY ANALYSIS');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end

    % Create output directories if non existent 
    if ~exist(path_data_out, 'dir'); mkdir(path_data_out); end
    if ~exist(path_img_out,'dir'); mkdir(path_img_out); end

    % Run script 
    s03_connectivity_reliability_analysis;
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
    
end

% -------------------------------------------------
% Correlation analysis  
% -------------------------------------------------

if flag.correlation_analysis
    
    % Define input/output data/paths
    data_out = filename.correlation;
    path_eeg_in = path.eeg_feature;
    path_bold_in = path.bold_preproc;
    path_data_out = path.correlation;
    path_img_out = strcat('IMAGES/', path_data_out);
    path_report_out = path.report;

    % Create report object
    % and open report 
    if flag.report ~= 0
        
        my_report = 'RESULTS - CORRELATION ANALYSIS';
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end
    
    % Go through sessions 
    for se = 1 : length(sessions)
        
        % Session 
        session = sessions(se);
        
        % Run script
        s02_correlation_analysis;
    
    end % sessions 
    
    % Close report 
    if flag.report ~= 0
        close(R);
    end
            
end

% -------------------------------------------------
% Group correlation analysis (ALL)
% -------------------------------------------------

if flag.group_correlation_analysis

    % Define input/output data/paths 
    data_in = filename.correlation;
    data_out = filename.correlation_group;
    path_data_in = path.correlation;
    path_data_out = path.correlation_group;
    path_img_out =  strcat('IMAGES\', path_data_out);
    path_report_out = path.report;
    
    % Create report object
    % and open report 
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
    if ~exist(path_data_out, 'dir'); mkdir(path_data_out); end
    if ~exist(path_img_out,'dir'); mkdir(path_img_out); end

    % Run script 
    s03_group_correlation_analysis;
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
    
end

% -------------------------------------------------
% Group correlation analysis (SESSIONS) 
% -------------------------------------------------

if flag.group_correlation_analysis_ses 

    % Define input/output data/paths 
    data_in = filename.correlation;
    data_out = filename.correlation_group;
    path_data_in = path.correlation;
    path_data_out = path.correlation_group_ses;
    path_img_out =  strcat('IMAGES\', path_data_out);
    path_report_out = path.report;
    
    % Create report object
    % and open report 
    if flag.report ~= 0
        
        my_report = strcat('RESULTS - GROUP', ...
            ' CORRELATION ANALYSIS (SESSIONS)');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\', my_report);
        open(R)
    
    end

    % Run script 
    s03_group_correlation_analysis_ses;
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
    
end

% -------------------------------------------------
% Correlation reliability analysis 
% -------------------------------------------------

if flag.correlation_reliability_analysis

    % Define input/output data/paths 
    data_in = filename.correlation;
    data_in_power = filename.power;
    data_in_power_reliab = filename.power_reliability;
    data_in_connectivity = filename.connectivity;
    data_in_connectivity_reliab = filename.connectivity_reliability;
    data_out = filename.correlation_reliability;
    path_data_in = path.correlation;
    path_data_in_power = path.power;
    path_data_in_power_reliab = path.power_reliability;
    path_data_in_connectivity = path.connectivity; 
    path_data_in_connectivity_reliab = path.connectivity_reliability;
    path_data_out = path.correlation_reliability;
    path_img_out =  strcat('IMAGES\', path_data_out);
    path_report_out = path.report;
    
    % Create report object
    % and open report 
    if flag.report ~= 0
        
        my_report = strcat('RESULTS -', ...
            ' CORRELATION RELIABILITY ANALYSIS');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
    
    end

    % Create output directories if non existent 
    if ~exist(path_data_out, 'dir'); mkdir(path_data_out); end
    if ~exist(path_img_out,'dir'); mkdir(path_img_out); end

    % Run script 
    s03_correlation_reliability_analysis;
    
    % Close report 
    if flag.report ~=0
        close(R);
    end
    
end

% -------------------------------------------------
% Estimate order of the ACF model 
% -------------------------------------------------
if flag.estimate_acf_order
    
    % If cross-validation is to be performed
    % across sessions, it won't be necessary 
    % to compute the order of the ACF
    if ~strcmp(cv_method, 'sessions')
        
        % Define input/output data/paths 
        data_in = filename.bold_preproc;
        data_deconv_in = filename.bold_deconv;
        path_data_in = path.bold_preproc;
        path_deconv_in = path.bold_deconv;
        path_img_out = strcat('\IMAGES\', ...
            path_pars);

       % Go through sessions 
        for se = 1 : length(sessions)
        
            % Session 
            session = sessions(se);
        
            % Run script 
            s04_estimate_acf_order;
            
        end % sessions 
        
    end
    
end

% -------------------------------------------------
% Optimize Cross-Validation Parameters 
% -------------------------------------------------
if flag.optimize_cv_pars
   
    % If cross-validation is to be performed
    % across sessions, it won't be necessary 
    % to compute the optimal cv parameters 
    if ~strcmp(cv_method, 'sessions')
        
        % Define input/output data/paths 
        path_eeg_in = path.eeg_feature;
        path_bold_in = path.bold_preproc; 

        % Optimize CV parameters 
        % for current CV method 
        s04_optimize_cv_pars;
        
    end
        
end
    
% -------------------------------------------------
% Fit EEG-BOLD models
% -------------------------------------------------
if flag.fit_models

    % Define input/output data/paths 
    path_eeg_in = path.eeg_feature;
    path_bold_in = path.bold_preproc;
    path_data_out = path.model; 
    path_img_out = strcat('IMAGES\', path_data_out);
    path_report_out = path.report;

    if strcmp(cv_method, 'sessions')
        
        % Run script
        s05_fit_models
    
    else
        
        % Go through sessions if CV is 
        % not to run across sessions 
        for se = 1 : length(sessions)

            session = sessions(se);

            % Fit model for current 
            % CV method 
            s05_fit_models;

        end % sessions
        
    end 
    
end

% -------------------------------------------------
% Report EEG-BOLD Models
% -------------------------------------------------
if flag.report_models
    
    % Define input/output data/paths  
    path_eeg_in = path.eeg_feature;
    path_bold_in = path.bold_preproc;
    path_data_in = path.model;
    
    path_img_out = strcat('IMAGES\', path_data_in);
    path_report_out = path.report;
    
    % Create the report object
    my_report = strcat('RESULTS -', ...
        " ", upper(reg_models), ' MODELS');
    R = Report(my_report, 'pdf');
    R.Layout.Landscape = true;
    R.OutputPath = strcat(path_report_out, ...
        '\', my_report);
    open(R)
    
    % Go through sessions if CV is 
    % not to run across sessions 
    for se = 1 : length(sessions)

        session = sessions(se);

        % Run script 
        s05_report_models;

    end % sessions 

    close(R)
    
end

% CODE NEEDS TO BE CHANGED FROM HERE TO ACCOMODATE RUNS 

% -------------------------------------------------
% Group Statistics of EEG-BOLD Models 
% -------------------------------------------------

if flag.group_model_stats
    
    % Define input/output data/paths  
    data_out = filename.model_group;
    path_data_in = path.model;
    path_data_out = path.model_group; 
    
    path_img_out = strcat('IMAGES\',path_data_out);
    path_report_out = path.report;
    
    % Go through CV methods 
    for r = 1 : length(reg_models)
        
        reg_model = reg_models(r);

        if flag.report ~= 0
            
            % Create the report object 
            % and open report 
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
        
        % Close report 
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
        % and open report 
        my_report = strcat('RESULTS - ', ...
            ' MODELS PERFORMANCE COMPARISON');
        R = Report(my_report,'pdf');
        R.Layout.Landscape = true;
        R.OutputPath = strcat(path_report_out, ...
            '\',my_report);
        open(R)
  
     end
     
    % Run script 
    s07_compare_model_performance;
    
    % Close report 
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
