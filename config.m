% Define configuration parameters for the analysis

%------------------------------------------------------------------
% Pipeline Parameters 
%------------------------------------------------------------------

% If a variable here defined is a list, then all the scripts will 
% loop through that list. If a variable here defined is a string, 
% then this is a fixed parameter for the entire pipeline. 

% Subjects 
subjects = ["sub-patient002", "sub-patient003", ...
       "sub-patient005", "sub-patient006", "sub-patient007", ...
       "sub-patient008", "sub-pilot011", "sub-pilot015", ...
       "sub-pilot018"];

% Task          -
% 'task-rest','task-calib'
task = "task-rest"; 

% Sub-task
% 'sub_task-EO','sub_task-EC'
sub_task = '';           

% BOLD RSN extraction method
% 'ic_dmn','avg_dmn'
rsn_method = "ic_dmn";              

% EEG TF-decomposition method
% 'wavelet','welch'
tf_method = "wavelet";               

% BOLD deconvolution method
% 'voxel_wise','time_series'
deconv_method = "time_series";      

% EEG feature decomposition metrics
% 'lc4','lc6','rmsf','tp'
metrics = ["lc4", "lc6", "rmsf", "tp"];

% Regression models 
% 'l21_1','elasticnet'
reg_models = "l21_1";   

% Cross-validation method
% 'nondep','regular','blocked'
cv_method = "nondep";     

% Surrogate method for statistical
% filtering of the connectivity
% 'block_shuffle', 'phase_shuffle'
surrogate_method = 'block_shuffle';

% Threshold of the DMN mask
dmn_thr = 1;

% Default channel and
% delay for plotting 
plotting_channel = 10;
plotting_delay = 3;

%------------------------------------------------------------------
% Analysis Parameters 
%------------------------------------------------------------------

% These are the parameters that either:
% 1) Are used in more than one script
% 2) Need to be reported in the Methods section of the paper

% Sampling frequency (Hz) 
fs_eeg = 250;               % EEG original sampling frequency
fs_bold = 1/1.26;           % BOLD original sampling frequency
fs_analysis = 4;            % Analysis intermediate sampling frequency

% Frequency range (Hz)
f_min = 1;
f_max = 30;
n_freq = 100;

% EEG filters (Hz)
highpass_filter = 1;          
lowpass_filter = 40;         

% Windows for TF decomposition 
tf_sliding_window_seconds = 4;  
tf_wavelet_kernel_seconds = 2; 

% Number of windows for 
% Welch's method 
n_wins_welch = 8;

% Number of surrogates for 
% statistical filtering of
% the connectivity estimates 
n_surrogates = 300;

% Window for HRF convolution 
hrf_kernel_seconds = 32;       

% Supported power and connectivity metrics 
power_metrics = ["lc4","lc6","rmsf","tp"];
connectivity_metrics = ["icoh_wnd", "icoh_cc"];

% Confidence level for the
% auto-regressive model used
% to model the EEG and BOLD
% time-series 
acf_conf = 0.99;

% Dimensions of BOLD 
% images 
n_xvoxs = 100;
n_yvoxs = 100;
n_zvoxs = 60;

% Group correlation 
% tstat threshold 
thresh_corr = 1;

% Group models
% tstat threshold
thresh_model = 1.5; 

% Number of global
% randomization tests to do 
n_rand = 1000;

%------------------------------------------------------------------
% EEG Events 
%------------------------------------------------------------------

% Define event markers for each task 
switch task
    
    case 'task-rest'
        
        markers_task = 'Scan Start';
        
    case 'task-calib'
        
        markers_task = 'Scan Start';
        markers_task_start = 'Scan Start';
        markers_task_stop = 'LR'; 
        
        switch sub_task
            
            case 'sub_task-EO'
                
                markers_sub_task_start = ["EO", "Rest"];
                markers_sub_task_stop = "EC"; 
                sub_task_order = 1;
                
            case 'sub_task-EC'
                
                markers_sub_task_start = "EC";
                markers_sub_task_stop = ["EO", "Rest"];
                sub_task_order = 2;
                
        end
        
end

%------------------------------------------------------------------
% Filename templates 
%------------------------------------------------------------------

% Rules for naming files - 1) modality
%                          2) processing underwent                          

% EEG/BOLD Raw data
filename.eeg_raw =              "ICremoved.set";       
filename.bold_raw =             strcat(rsn_method,'.txt');
filename.bold_img_raw =         "filtered_func_data_preprocessed.nii.gz";

% EEG/BOLD Markers
filename.eeg_markers =          strcat(task,'_timing_file.txt');
filename.eeg_markers_sub_task = strcat(sub_task,'_timing_file.txt');
filename.bold_markers =         strcat(task,'_timing_file.txt');
    
% EEG/BOLD Processed
filename.eeg_processed =        "eeg_processed.mat";
filename.bold_processed =       "bold_processed.txt";
filename.bold_img_processed =   "bold_img_processed.txt";

% EEG/BOLD Derivatives
filename.eeg_feature =          "eeg_feature.txt";
filename.eeg_feature_eeg_fs =   "eeg_feature_eeg_fs.txt";
filename.eeg_feature_conv =     "eeg_feature_conv.txt";
filename.eeg_feature_delay =    "eeg_feature_delay.txt";
filename.bold_deconv =          "bold_processed_deconv.txt";
filename.bold_img_deconv =      "bold_img_processed_deconv.txt";

% BOLD masks of the DMN 
filename.bold_mask_dmn =        "mask_dmn.nii.gz";
filename.bold_mask_dmn_bin =    "mask_dmn.txt";

% EEG Frequency  
filename.frequency =            "";

% EEG-BOLD Correlation
filename.correlation =          "corr_stats.mat";
filename.correlation_group =    "corr_gstats.mat";

% EEG-BOLD Model
filename.model_group =          "model_gstats.mat";

% EEG-BOLD Model 
% Performance Comparisson 
filename.anova =                "anova.mat";

%------------------------------------------------------------------
% Filename templates 
%------------------------------------------------------------------

% Directories rules - 1) main path, followed by:
%                           1.1) DATA           minimal processing
%                           1.2) DERIVATIVES    altered structure 
%                           1.3) RESULTS        draw conclusions from
%                           1.4) IMAGES         images  
%                           1.5) REPORTS        automatic reports 
%                     2) followed by the dataset path
%                     3) followed by a specific path 

% Define the main path 
path.main = 'C:\Users\marta\Documents\LASEEB\MigN2Treaty';

% Define path that defines the dataset
% (Common to almost every path) 
data_path = {subjects,task,sub_task};

% Define path that defines the set of methods specified 
method_path = {rsn_method,tf_method,reg_models,cv_method};

% EEG/BOLD Raw data
path.eeg_raw =              strcat('DATA\',fullfile(data_path{:}),'\eeg');
path.bold_raw =             strcat('DATA\',fullfile(data_path{:}),'\func');
path.bold_img_raw =         strcat('DATA\',fullfile(data_path{:}),'\func');

% EEG/BOLD markers
path.eeg_markers =          strcat('DATA\',fullfile(data_path{:}),'\eeg');
path.bold_markers =         strcat('DATA\GROUP\',fullfile...
                            (data_path{2:3}),'\eeg');

% EEG/BOLD Processed 
path.eeg_processed =        strcat('DATA\',fullfile(data_path{:}),'\eeg');
path.bold_processed =       strcat('DATA\',fullfile(data_path{:}),'\func\', ...
                            rsn_method);
path.bold_img_processed =   strcat('DATA\',fullfile(data_path{:}),'\func');                        

% EEG/BOLD Derivatives 
path.eeg_feature =          strcat('DERIVATIVES\',fullfile(data_path{:}), ... 
                            '\eeg\',method_path{2});
path.bold_deconv =          strcat('DERIVATIVES\',fullfile(data_path{:}), ...
                            '\func\',method_path{1});
path.bold_img_deconv =      strcat('DERIVATIVES\',fullfile(data_path{:}), ...
                            '\func\',method_path{1});

% EEG Frequency 
path.frequency =            strcat('DERIVATIVES\',fullfile(data_path{:}), ...
                            '\eeg\',method_path{2},'\frequency_analysis\');

% EEG-BOLD Correlation 
path.correlation =          strcat('RESULTS\',fullfile(data_path{:}), ...
                            '\correlation_analysis\', ...
                        	fullfile(method_path{1:2}));
path.correlation_group =   strcat('RESULTS\GROUP\', ...
                            fullfile(data_path{2:3}), ...
                            '\correlation_analysis\', ...
                            fullfile(method_path{1:2}));

% EEG-BOLD Model
pre =                       repmat(strcat('RESULTS\', ...
                            fullfile(data_path{:})', ...
                            '\models\'),1,length(reg_models));
suf =                       repmat(fullfile(method_path{1:3}), ...
                            length(subjects),1);
path.model =                strcat(pre,suf);
path.model_group =          strcat('RESULTS\GROUP\', ...
                            fullfile(data_path{2:3}), ...
                            '\models\',fullfile(method_path{:}));
                        
                        
% Model Performance 
cat =                       [fullfile(method_path{:}) strcat(fullfile...
                            (method_path{1:2}),'\',strjoin(reg_models, ...
                            '_vs_'),'\',cv_method)];
path.compare_performance =  strcat('RESULTS\GROUP\', ...
                            fullfile(data_path{2:3}), ...
                            '\performance\',cat);                                       
                        
% Report
path.report =               strcat('REPORTS\',fullfile(data_path{2:3}), ...
                            '\',fullfile(method_path{1:2}));
                        
% Parameters                         
path.pars =                 strcat('PARS\',fullfile(data_path{2:3}), ...
                            '\',fullfile(method_path{1:2}));

%------------------------------------------------------------------
% Image settings 
%------------------------------------------------------------------

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','off');
