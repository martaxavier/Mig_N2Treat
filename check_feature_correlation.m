% Checks correlation among imput features 

clear all

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','off');

% Subjects and task
subjects = ["sub-patient002", "sub-patient003", ...
     "sub-patient005", "sub-patient006", "sub-patient007", ...
     "sub-patient008"];
 task = 'task-rest';

% Define paths
path.main = 'C:\Users\marta\Documents\LASEEB\MigN2Treaty';
path_X = fullfile('DERIVATIVES', subjects, task, 'eeg\wavelet', 'analytical');
path_Y = fullfile('DERIVATIVES', subjects, task, 'eeg\wavelet', 'analytical');
path_Z = fullfile('DERIVATIVES', subjects, task, 'eeg\wavelet', 'analytical');
path_data_out = fullfile('DERIVATIVES', subjects, task, 'eeg\wavelet', 'compare_features');
path_img_out = fullfile('IMAGES',path_data_out);

% Define data
data_in_X = 'icoh_wnd_eeg_feature.txt';
data_in_Y = 'icoh_wne_eeg_feature.txt';
data_in_Z = 'icoh_bc_eeg_feature.txt';
id_X = 'NodeDegree';
id_Y = 'NodeEfficiency';
id_Z = 'BetweennessC';
%id_Y = 'SurrogatePhaseShuffle';

% Get metric parameters
metric = 'icoh';
get_metric_pars;
n_subjects = length(subjects);

 for s = 1 : n_subjects
    
    % Create output directories if non-existent 
    if ~exist(path_img_out(s), 'dir'); mkdir(path_img_out(s)); end

    con_X = dlmread(fullfile(path_X(s), data_in_X));
    con_Y = dlmread(fullfile(path_Y(s), data_in_Y));
    con_Z = dlmread(fullfile(path_Z(s), data_in_Z));
    con_X = reshape(con_X, [size(con_X, 1) n_chans n_bands]);
    con_Y = reshape(con_Y, [size(con_Y, 1) n_chans n_bands]);
    con_Z = reshape(con_Z, [size(con_Z, 1) n_chans n_bands]);
    
    for c = 1 : n_chans
        
        fig = figure('Name', 'Correlation Matrices');
        
        for b = 3 : 3
            %subplot(2, 2, b);
            XY = [con_X(:, c, b) con_Y(:, c, b) con_Z(:, c, b)];
            corrplot(XY, 'varNames', {id_X, id_Y id_Z});
            th = findall(gcf, 'type', 'text', 'String', '{\bf Correlation Matrix}'); 
            th.String = strcat('CorrCoef,', " " , id_bands(b), ' band');
        end
        
        img_out = strcat(id_X, '_vs_', id_Y, '_vs_', id_Z, '_chan', ...
            num2str(c), '_', id_bands(b), '.png');
        saveas(fig,fullfile(path_img_out(s), img_out));
        
    end
    
 end