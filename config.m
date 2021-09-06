% Define configuration parameters for the analysis

%------------------------------------------------------------------
% Pipeline Parameters 
%------------------------------------------------------------------

% If a variable here defined is a list, then all the scripts will 
% loop through that list. If a variable here defined is a string, 
% then this is a fixed parameter for th0e entire pipeline.
% Each combination of all the variables in this section will originate
% results that are to be saved in a specific directory 

% Subjects 
%   subjects = ["sub-patient002", "sub-patient003", ...
%   "sub-patient005", "sub-patient006", "sub-patient007", ...
%   "sub-patient008", "sub-patient012"];
%   subjects = ["sub-32" "sub-35" "sub-36" "sub-37" "sub-38" ...
%      "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-45", ...
%       "sub-46" "sub-47" "sub-48" "sub-49" "sub-50"]; 
       subjects = ["sub-04" "sub-05" "sub-06" ...
           "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" ...
           "sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-19" ...
           "sub-20" "sub-21" "sub-22" "sub-23" "sub-24" "sub-25" ...
           "sub-27" "sub-29"];
      
% Dataset
% 'Mig_N2Treat', 'NODDI', 'PARIS';
dataset = 'PARIS'; 

% Task   
% 'task-rest','task-calib'
task = "task-rest"; 

% Sub-task
% 'sub_task-EO','sub_task-EC'
sub_task = '';    

% Sessions
% "run-1", "run-2", "run-3"
sessions = ["run-2"];

% BOLD RSN extraction method
% 'ic_dmn','avg_dmn'
rsn_method = "ic_dmn_group";              

% EEG TF-decomposition method 
% 'wavelet','welch'
tf_method = "wavelet";    

% Method for statistical filtering 
% of the connectivity matrices 
% 'analytical', 'surrogate', 'none' 
stat_filt_method = 'analytical';

% Method for generation of the 
% connectivity surrogates 
% 'block_shift', 'phase_shuffle', ''
surr_method = ''; 

% BOLD deconvolution method
% 'voxel_wise','time_series'
deconv_method = "time_series";      

% EEG feature decomposition metrics
% 'lc4','lc6','rmsf','tp','icoh_wnd','icoh_bc','lc5'
metrics = ["lc4"];

% Reliability metric
reliab_metric = "icc";

% Regression models 
% 'l21_1', 'l2_1'
reg_models = "l2_1";   

% Hierarchical clustering 
% distance metric and algorithm
hc_distance = ["average", "correlation"];
    
% Cross-validation method
% 'nondep','regular','blocked','sessions','one_class'
cv_method = "blocked";     

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
fs_bold = 1/2;              % BOLD original sampling frequency
fs_analysis = 4;            % Analysis intermediate sampling frequency

% Frequency range (Hz)
f_min = 1;
f_max = 30;
n_freq = 30;

% EEG filters (Hz)
highpass_filter = 1;          
lowpass_filter = 40;         

% Windows for TF decomposition 
tf_sliding_win_seconds = 4;  
tf_wavelet_kernel_seconds = 2; 

% Number of windows for 
% Welch's method 
n_wins_welch = 8;

% Number of surrogates for 
% statistical filtering of
% the connectivity estimates 
n_surrs = 100;

% Window for HRF convolution 
hrf_kernel_seconds = 32;       

% Supported power and connectivity metrics 
power_metrics = ["lc4" "lc5" "lc6" "rmsf" "tp"];
connectivity_metrics = ["icoh_wnd" "wpli_wnd", ...
    "icoh_wne" "wpli_wne" "icoh_bc" "wpli_bc"];

% Confidence level for the
% auto-regressive model used
% to model the EEG and BOLD
% time-series 
acf_conf = 0.95;

% Dimensions of BOLD 
% images 
n_xvoxs = 100;
n_yvoxs = 100;
n_zvoxs = 60;

% Group model/correlation 
% tstat threshold (plot)
thresh_tstat = 1.5;

% False discovery rate 
% threshold 
thresh_fdr = 0.05;

% Number of randomnization tests
% for topographic consistency 
n_rand = 1000;

% Number of iterations for split 
% half reliability test
n_iter_reliab = 100;

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
filename.eeg_preproc =        "eeg_preproc.mat";
filename.bold_preproc =       "bold_preproc.txt";
filename.bold_img_preproc =   "bold_img_preproc.txt";

% EEG/BOLD Derivatives
filename.eeg_feature =          "eeg_feature.txt";
filename.eeg_feature_eeg_fs =   "eeg_feature_eeg_fs.txt";
filename.eeg_feature_conv =     "eeg_feature_conv.txt";
filename.eeg_feature_delay =    "eeg_feature_delay.txt";
filename.bold_deconv =          "bold_preproc_deconv.txt";
filename.bold_img_deconv =      "bold_img_preproc_deconv.txt";

% BOLD masks of the DMN 
filename.bold_mask_dmn =        "mask_dmn.nii.gz";
filename.bold_mask_dmn_bin =    "mask_dmn.txt";

% EEG Power  
filename.power =                "power.mat";
filename.power_reliability =    "power_reliab.mat"; 

% EEG Static Connectivity 
filename.connectivity =             "connectivity.mat";
filename.connectivity_reliability = "connectivity_reliability.mat";

% EEG-BOLD Correlation
filename.correlation =              "corr_stats.mat";
filename.correlation_group =        "corr_gstats.mat";
filename.correlation_reliability =  "corr_reliab.mat";

% EEG-BOLD Model
filename.model_group =          "model_gstats.mat";

% EEG-BOLD Model 
% Performance Comparisson 
filename.anova =                "anova.mat";

%------------------------------------------------------------------
% Path templates 
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
data_path = {subjects, task, sub_task, sessions};

% Define path that defines the set of methods specified 
%stat_method = fullfile(stat_filt_method, surr_method);
method_path = {rsn_method, tf_method, reg_models, cv_method};

% EEG/BOLD Raw data
pre =                       repmat(fullfile(data_path{1:3}), max(length(data_path{4}), 1), 1)';
suf =                       repmat(data_path{4}, max(length(data_path{1}), 1), 1);
path.eeg_raw =              strcat(dataset, '\DATA\eeg\', fullfile(pre, suf));
path.bold_raw =             strcat(dataset, '\DATA\func\', fullfile(pre, suf));
path.bold_img_raw =         strcat(dataset, '\DATA\func\', fullfile(pre, suf));

% EEG/BOLD markers
path.eeg_markers =          strcat(dataset, '\DATA\eeg\', fullfile(pre, suf));
path.bold_markers =         strcat(dataset, '\DATA\eeg\', fullfile(data_path{2:3}), '\GROUP');

% EEG/BOLD Processed 
path.eeg_preproc =        strcat(dataset, '\DATA\eeg\', fullfile(pre, suf));
path.bold_preproc =       strcat(dataset, '\DATA\func\', fullfile(pre, suf, rsn_method));
path.bold_img_preproc =   strcat('DATA\func\', strcat(pre, suf));                        

% EEG/BOLD Derivatives 
path.eeg_feature =          strcat(dataset, '\DERIVATIVES\eeg\', fullfile(pre, suf, '\features\', fullfile(method_path{2})));
path.bold_deconv =          strcat(dataset, '\DERIVATIVES\func\', fullfile(pre, suf, method_path{1}));
path.bold_img_deconv =      strcat(dataset, '\DERIVATIVES\func\', fullfile(pre, suf, method_path{1}));

% EEG Power  
path.power =                fullfile(dataset, '\DERIVATIVES\eeg\', fullfile(pre, suf), '\power_analysis\');
path.power_reliability =    fullfile(dataset, '\DERIVATIVES\eeg\GROUP\', fullfile(data_path{2:3}), '\power_reliability_analysis\', upper(reliab_metric));

% EEG Static Connectivity 
path.connectivity =             fullfile(dataset, '\DERIVATIVES\eeg\', fullfile(pre, suf), '\connectivity_analysis\');
path.connectivity_reliability = fullfile(dataset, '\DERIVATIVES\eeg\GROUP\', fullfile(data_path{2:3}), '\connectivity_reliability_analysis\', upper(reliab_metric));

% EEG-BOLD Correlation 
path.correlation =              strcat(dataset, '\RESULTS\', fullfile(pre, suf), '\correlation_analysis\', fullfile(method_path{1:2}));
path.correlation_group =        strcat(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:3}), '\correlation_analysis\', fullfile(method_path{1:2}));
path.correlation_group_ses =    strcat(dataset, '\RESULTS\', fullfile(data_path{1:3}), '\correlation_analysis_ses\', fullfile(method_path{1:2}));    
path.correlation_reliability =   fullfile(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:3}), '\correlation_reliability_analysis\', fullfile(method_path{1:2}), upper(reliab_metric));

% EEG-BOLD Model
path.model =                    strcat(dataset, '\RESULTS\', fullfile(pre, suf), '\models\', fullfile(method_path{1:3}));
path.model_group =              strcat(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:4}), '\models\', fullfile(method_path{1:3}));                        

path.model_ses =                strcat(dataset, '\RESULTS\', fullfile(data_path{1:3}), '\models\', fullfile(method_path{1:3}));
path.model_ses_group =          strcat(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:3}), '\models\', fullfile(method_path{1:3}));                        

path.model_hc =                 strcat(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:3}), '\hierarchical_clustering\', fullfile(method_path{1:3}));
path.model_one_class =          strcat(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:3}), '\one_class_models\', fullfile(method_path{1:3}));

% Model Performance 
ccat =                       [fullfile(method_path{:}), strcat(fullfile(method_path{2:3}), '\', strjoin(reg_models, '_vs_'), '\', cv_method)];
path.compare_performance =  strcat(dataset, '\RESULTS\GROUP\', fullfile(data_path{2:3}), '\performance\', ccat);                                       
                        
% Report
path.report =               strcat(dataset, '\REPORTS\', fullfile(data_path{2:3}), '\', fullfile(method_path{1:2}));
                        
% Parameters                         
path.pars =                 strcat(dataset, '\PARS\', fullfile(data_path{2:3}), '\', fullfile(method_path{1:2}));

%------------------------------------------------------------------
% Image settings 
%------------------------------------------------------------------

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','off');
