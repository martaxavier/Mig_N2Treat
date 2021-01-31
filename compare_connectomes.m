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
path_Y = fullfile('DERIVATIVES', subjects, task, 'eeg\wavelet', 'surrogate\phase_shuffle');
path_data_out = fullfile('DERIVATIVES', subjects, task, 'eeg\wavelet', 'compare_connectomes');
path_img_out = fullfile('IMAGES',path_data_out);

% Define data
data_in_X = 'conspec_topo_filtered_icoh.mat';
data_in_Y = 'conspec_topo_filtered_icoh.mat';
id_X = 'Analytical'; id_Y = 'SurrogatePhaseShuffle';

% Get parameters
metric = 'icoh';
get_metric_pars;

n_subjects = length(subjects);
sim_corr = zeros(n_subjects, 1);
sim_dice = sim_corr;

 for s = 1 : n_subjects

    % Create output directories if non-existent 
    if ~exist(path_img_out(s), 'dir'); mkdir(path_img_out(s)); end
    
    load(fullfile(path_X(s), data_in_X)); con_X = conspec_topo;
    load(fullfile(path_Y(s), data_in_Y)); con_Y = conspec_topo; 
    sim_corr(s) = connectome_similarity(con_X, con_Y, 'corr_coef');
    sim_dice(s) = connectome_similarity(con_X, con_Y, 'dice_coef');
    
    % Plot a random time-point and band, and the 
    % Dice and Pearson coefficient for current subject 
    
    for t = 1 : 2
        
        if t == 1
            rand_pnt = 100; 
        else
            rand_pnt = 1000;
        end
    
        for b = 1 : n_bands 

            my_title = strcat('Connectome at time-point', " ", ...
                num2str(rand_pnt),' , for the', " ", ...
                id_bands(b), ' band');
            fig = figure('Name', my_title);
            fig.Position(3 : 4) = fig.Position(3 : 4)*5;

            subplot(1, 2, 1);
            imagesc(squeeze(con_X(rand_pnt, :, :, b))); colorbar; 
            title(strcat(id_X, ' (corrCoef =', " ", num2str(sim_corr(s)), ...
                '), (diceCoeff =', " ", num2str(sim_dice(s)),')'), 'FontSize', 25);
            xlabel('Channels','FontSize', 24); 
            ylabel('Channels','FontSize', 25);
            xticks(1 : n_chans); yticks(1 : n_chans);
            xticklabels(cellstr(id_chans));
            yticklabels(cellstr(id_chans));

            subplot(1, 2, 2);
            imagesc(squeeze(con_Y(rand_pnt, :, :, b))); colorbar; 
            title(strcat(id_Y, ' (corrCoef =', " ", num2str(sim_corr(s)), ...
                '), (diceCoeff =', " ", num2str(sim_dice(s)),')'), 'FontSize', 25);
            xlabel('Channels','FontSize', 24); 
            ylabel('Channels','FontSize', 25);
            xticks(1 : n_chans); yticks(1 : n_chans);
            xticklabels(cellstr(id_chans));
            yticklabels(cellstr(id_chans));

            set(gcf, 'PaperUnits', 'inches');
            x_width = 30; y_width = 11;
            set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

            % Save generated images in output path 
            img_out = strcat(upper(metric), '_', id_X, ...
                '_vs_', id_Y, '_Time', num2str(rand_pnt), ...
                '_', id_bands(b), '.png');
            saveas(gcf,fullfile(path_img_out(s), img_out));
            
        end % finish looping through bands
     
    end % finish looping through time samples 
    
 end % finish looping through subhects 
 
 % Cross-subject similarity, X
 % Cross-subject similarity, Y

    function [sim] = connectome_similarity(con_X, con_Y, sim_measure)

    % Compares the connectivity matrices 'con_X' and 'con_Y' using a 
    % similarity measure specified by 'sim'
    % 
    % INPUTS
    %   con_X, con_Y        input connectivity matrices
    %   sim                 similarity measure
    %                           'corr_coef' - Pearson's correlation 
    %                           'dice_coef' - Dice coefficient
    %                           'manhattan' - Manhattan (L1) distance
    %                           'euclidean' - Euclidean (L2) distance 
    %                           the latter two are applied onto the 

    switch sim_measure
        % Pearson's correlation between
        % each entry of the connecivity matrix
        case 'corr_coef'         
            sim = corrcoef(con_X(:), con_Y(:));
            sim = sim(1, 2);
        % Dice coefficinent of the entire 
        % dynamic connectome 
        case 'dice_coef'
            sim = dice(logical(con_X), logical(con_Y));
    end

    end