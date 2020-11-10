close all
path='C:\Users\marta\Documents\LASEEB\MigN2Treat'; cd(path)

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','off');

% ------------------------------------------------------------
% Input configuration parameters 
% ------------------------------------------------------------              

% Write configuration struct
cfg.shift = 'conv';                 % 'conv','delay'
cfg.cv = 'nondep';                  % 'regular','nondep','blocked'

cfg.metric = 'lc';                  % 'lc','rmsf','tp'
cfg.regress = 'l21_1';              % 'l21_1','elasticnet'
cfg.compare = 'regress';             % 'metric','regress'
cfg.stats = 'anova';              % 'anova','profile'
cfg.write = 1;                      % 1 to write, 0 to don't

cfg.subjects = ["sub-patient002","sub-patient003",...
    "sub-patient006","sub-pilot011","sub-pilot015",...
    "sub-pilot018"];

% Load the .mat file containing
% model configuration parameters 
% for each of the EEG-BOLD models used 
load('modelpars','modelpars');

% ------------------------------------------------------------
% Write and/or read group stats configuration file 
% ------------------------------------------------------------

% Write group stats file with model labels and paths 
if cfg.write; write_group_stats_file(cfg); end

% Read group stats file with the labels and
% paths of models to compare 
switch cfg.stats
    
    case 'profile'
        
        paths_in = strcat(cfg.regress,'/.mat/run_',...
            cfg.metric,'_',cfg.shift,'_',cfg.cv,'.mat');
        
    case 'anova'
        
        file_id = fopen(strcat('RESULTS/',...
            'group_stats/group_stats.txt'),'r');
        config = textscan(file_id,'%s %s\n'); 
        fclose(file_id);
        labels_in = string(config{1,1}); 
        paths_in = string(config{1,2});
        
end

switch cfg.stats
    case 'profile'
        
        % Create and/or specify output directory 
        path_out = strcat(path,'\RESULTS\',...
            '\group_stats\model_stats\',cfg.stats,'\',...
            cfg.regress,'\',cfg.cv,cfg.metric);
        
        % Read model parameters for current metric
        pars = modelpars.(upper(cfg.metric));
    
    case 'anova'
        
        % Create and/or specify output directory 
        path_out = strcat(path,'\RESULTS\',...
            '\group_stats\model_stats\',cfg.stats,'\',cfg.compare);
        
end 

if ~exist(path_out, 'dir'); mkdir(path_out); end


% ------------------------------------------------------------
% Go through subjects and 
% ------------------------------------------------------------

clear optimals; clear models;

% Go through subjects to read models
for s = 1 : length(cfg.subjects)
    
    subject = cfg.subjects(s);
    
    % Load models for current subject 
    path_run_in = strcat(path,'\RESULTS\',...
        subject,'\models\',paths_in);
    
    for r = 1 : size(path_run_in,1)
        load(path_run_in(r),'model','optimal')
        optimals(r,s) = optimal;
        models(r,s) = model;
    end
    
end

% ------------------------------------------------------------
% Perform anova or profile analysis
% ------------------------------------------------------------

switch cfg.stats
    case 'profile'
        get_profiles(models,pars,path,path_out);
    case 'anova'
        perform_anova(optimals,models,labels_in,path_out);
end

%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% perform_anova()             
% ============================================================

function [] = perform_anova(optimals,models,labels_in,path_out)

%   [] = perform_anova performs anova on the input data  
%        and plots and saves the corresponding boxplots 
%
%   INPUTS:
%
%     optimals         A struct containing the model results 
%                      at each cv fold for each subject and
%                      for each of the conditions being compared 
%     models           A struct containing the final model  
%                      results for each subject and for each 
%                      of the conditions being compared 
%   

% -------------------------------------------------
% Read input information  
% -------------------------------------------------

% Read size of the problem 
n_compare = size(optimals,1);
n_subjects = size(optimals,2);
n_folds = size(optimals(1,1).bic,1);

% -------------------------------------------------
% Write parameter matrices 
% -------------------------------------------------

% Pre-allocate parameters to evaluate 
bic = zeros(n_subjects*n_folds,n_compare);
nmse = bic; corr = bic;
df = zeros(n_subjects,n_compare);
lambda = df; rho = df; time = df;

for c = 1 : n_compare
    
    for s = 1 : n_subjects
        
        idxs = s*n_folds-n_folds+1:s*n_folds;
        bic(idxs,c) = optimals(c,s).bic;
        nmse(idxs,c) = optimals(c,s).nmse;
        corr(idxs,c) = optimals(c,s).corr;
        
        df(s,c) = models(c,s).df;
        lambda(s,c) = models(c,s).lambda;
        rho(s,c) = models(c,s).rho;
        time(s,c) = models(c,s).time;
        
    end
end

% -------------------------------------------------
% Perform anova  
% -------------------------------------------------

% bic 
[~,~,stats] = anova1(bic); c = multcompare(stats);
img_out = strcat('mulcmp_bic.mat');
save(fullfile(path_out,img_out),'c');

% nmse
[~,~,stats] = anova1(nmse); c = multcompare(stats);
img_out = strcat('mulcmp_nmse.mat');
save(fullfile(path_out,img_out),'c');

% corr
[~,~,stats] = anova1(corr); c = multcompare(stats);
img_out = strcat('mulcmp_corr.mat');
save(fullfile(path_out,img_out),'c');

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
saveas(gcf,fullfile(path_out,img_out),'fig');

% nmse
b = boxplot(nmse,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2});
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('NMSE','FontSize',22);
img_out = strcat('anova_nmse');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% corr
b = boxplot(corr,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2});
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('CORR','FontSize',22);
img_out = strcat('anova_corr');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% df
b = boxplot(df,'Colors',[0 0.4470 0.7410]);
set(b,{'linew'},{2});
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on';
set(gca,'XTickLabel',labels_in,'FontSize',18)
ylabel('DF','FontSize',22);
img_out = strcat('anova_df');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% -------------------------------------------------
% Plot results in table format!! 
% -------------------------------------------------
    
end

% ============================================================
% get_profiles()             
% ============================================================

function [] = get_profiles(models,pars,path,path_out)

%   [] = get_profiles performs analaysis of the frequency/
%        delay/channel profiles, plots and saves the results 
%
%   INPUTS:
%
%     optimals         A struct containing the model results 
%                      at each cv fold for each subject and
%                      for each of the conditions being compared 
%     models           A struct containing the final model  
%                      results for each subject and for each 
%                      of the conditions being compared 
%   

% -------------------------------------------------
% Read input information  
% -------------------------------------------------

% Save model parameters for current 
% metric as local variables 
save('pars.mat','-struct','pars');
load('pars.mat'); id_metric = id; 
clear pars id;

% Read size of the problem 
n_subjects = size(models,2);

% Pre-allocate efp matrix
efp = zeros(length(models(1).efp)-1,n_subjects);

% Go through subjects
for s = 1 : n_subjects
    efp(:,s) = models(s).efp(2:end);
end

efp_avg = mean(efp,2);
efp_abs = mean(abs(efp),2);

% Frequency labels
X_freq = categorical(cellstr(id_bands));
X_freq = reordercats(X_freq,cellstr(id_bands));

% Delay labels
X_del = categorical...
    (cellstr(id_delays));
X_del = reordercats...
    (X_del,cellstr(id_delays));


% -------------------------------------------------
% Frequency profiles 
% -------------------------------------------------

% Compute absolute frequency profiles 
efp_abs_freq = mean(reshape(efp_abs,...
    [n_chans*n_delays,n_bands]));

% Compute average frequency profiles 
efp_avg_freq = mean(reshape(efp_avg,...
    [n_chans*n_delays,n_bands]));

% Plot absolute frequency profiles 
figure_name = strcat...
    ('Absolute frequency profile across all subjects');    
figure('Name',figure_name);
bar(X_freq,efp_abs_freq,'FaceColor','#77AC30','BarWidth',0.6);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
img_out = strcat('PROFILE_ABSEFP_Bands');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot average frequency profiles 
figure_name = strcat...
    ('Average frequency profile across all subjects');    
figure('Name',figure_name);
bar(X_freq,efp_avg_freq,'FaceColor','#77AC30','BarWidth',0.6);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
img_out = strcat('PROFILE_AVGEFP_Bands');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Boxplot of absolute frequency profiles
efp_freq = permute(reshape(efp,...
    [n_chans*n_delays,n_bands,n_subjects]),[1 3 2]);
efp_freq = reshape(efp_freq,[n_chans*n_delays*n_subjects,4]);
boxplot(abs(efp_freq),'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',cellstr(X_freq),'FontSize',18)
ylabel('Absolute EFP','FontSize',22);
img_out = strcat('BOXPLOT_ABSEFP_Bands');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Boxplot of frequency profiles
boxplot(efp_freq,'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',cellstr(X_freq),'FontSize',18)
ylabel('EFP','FontSize',22);
img_out = strcat('BOXPLOT_AVGEFP_Bands');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% -------------------------------------------------
% Frequency profiles - of each delay 
% -------------------------------------------------

% Compute absolute frequency profiles at each delay
efp_abs_freqbydel = ...
    squeeze(mean(reshape(efp_abs,...
    [n_chans,n_delays,n_bands])));

% Compute average frequency profiles at each delay
efp_avg_freqbydel = ...
    squeeze(mean(reshape(efp_avg,...
    [n_chans,n_delays,n_bands])));

% Plot absolute frequency profiles at each delay
figure_name = strcat...
    ('Absolute frequency profile across all subjects, for each delay');    
figure('Name',figure_name); 
b = bar(X_del,efp_abs_freqbydel,'FaceColor','flat','BarWidth',0.9);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
for k = 1:size(efp_avg_freqbydel,2);b(k).CData = k; end
set(b,{'DisplayName'},cellstr(X_freq)'); legend('FontSize',16);
img_out = strcat('PROFILE_ABSEFP_BandbyDelay');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot average frequency profiles at each delay
figure_name = strcat...
    ('Average frequency profile across all subjects, for each delay');    
figure('Name',figure_name); 
b = bar(X_del,efp_avg_freqbydel,'FaceColor','flat','BarWidth',0.9);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
for k = 1:size(efp_avg_freqbydel,2);b(k).CData = k; end
set(b,{'DisplayName'},cellstr(X_freq)'); legend('FontSize',16);
img_out = strcat('PROFILE_AVGEFP_BandbyDelay');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% -------------------------------------------------
% Delay profiles 
% -------------------------------------------------

% Compute absolute delay profiles 
efp_abs_del = reshape(efp_abs,[n_chans,n_delays,n_bands]);
efp_abs_del = permute(efp_abs_del,[2,1,3]);
efp_abs_del = mean(reshape(efp_abs_del,...
    [n_delays, n_chans*n_bands]),2);

% Compute average delay profiles 
efp_avg_del = reshape(efp_avg,[n_chans,n_delays,n_bands]);
efp_avg_del = permute(efp_avg_del,[2,1,3]);
efp_avg_del = mean(reshape(efp_avg_del,...
    [n_delays, n_chans*n_bands]),2);

% Plot absolute delay profiles 
figure_name = strcat...
    ('Absolute delay profile across all subjects');    
figure('Name',figure_name);
bar(X_del,efp_abs_del,'FaceColor','flat','BarWidth',0.6);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
img_out = strcat('PROFILE_ABSEFP_Delay');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot average delay profiles 
figure_name = strcat...
    ('Average delay profile across all subjects');    
figure('Name',figure_name);
bar(X_del,efp_avg_del,'FaceColor','flat','BarWidth',0.6);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
img_out = strcat('PROFILE_AVGEFP_Delay');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot std of channel profiles 
efp_chan = reshape(efp,[n_chans,...
    n_delays,n_bands,n_subjects]);
efp_chan = reshape(efp_chan,[n_chans,...
    n_delays*n_bands*n_subjects]);
efp_chan_std = std(efp_chan,0,2);
efp_chan_std_abs = std(abs(efp_chan),0,2);

% Boxplot of absolute frequency profiles
efp_del = permute(reshape(efp,[n_chans,...
    n_delays,n_bands,n_subjects]),[1 3 4 2]);
efp_del = reshape(efp_del,[n_chans*n_bands*...
    n_subjects,n_delays]);
boxplot(abs(efp_freq),'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',cellstr(X_del),'FontSize',18)
ylabel('Absolute EFP','FontSize',22);
img_out = strcat('BOXPLOT_ABSEFP_Delay');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Boxplot of frequency profiles
boxplot(efp_del,'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',cellstr(X_del),'FontSize',18)
ylabel('EFP','FontSize',22);
img_out = strcat('BOXPLOT_AVGEFP_Delay');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% -------------------------------------------------
% Delay profiles - of each frequency  
% -------------------------------------------------

% Compute absolute frequency profiles at each delay
efp_abs_delbyfreq = ...
    squeeze(mean(reshape(efp_abs,...
    [n_chans,n_delays,n_bands])))';

% Compute average frequency profiles at each delay
efp_avg_delbyfreq = ...
    squeeze(mean(reshape(efp_avg,...
    [n_chans,n_delays,n_bands])))';

% Plot absolute frequency profiles at each delay
figure_name = strcat...
    ('Absolute delay profile across all subjects, for each frequency band');    
figure('Name',figure_name); 
b = bar(X_freq,efp_abs_delbyfreq,'FaceColor','flat','BarWidth',0.9);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
for k = 1:size(efp_abs_delbyfreq,2);b(k).CData = k; end
set(b,{'DisplayName'},cellstr(X_del)'); legend('FontSize',16);
img_out = strcat('PROFILE_ABSEFP_DelaybyBand');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot average frequency profiles at each delay
figure_name = strcat...
    ('Average delay profile across all subjects, for each frequency');    
figure('Name',figure_name); 
b = bar(X_freq,efp_avg_delbyfreq,'FaceColor','flat','BarWidth',0.9);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
for k = 1:size(efp_avg_delbyfreq,2);b(k).CData = k; end
set(b,{'DisplayName'},cellstr(X_del)'); legend('FontSize',16);
img_out = strcat('PROFILE_AVGEFP_DelaybyBand');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% -------------------------------------------------
% Channel profiles 
% -------------------------------------------------

% Get channel locations 
load(strcat(path,'\DATA\sub-patient002\',...
    'eeg\eeg_preprocessed.mat'),'EEG');
chanlocs = struct2cell(EEG.chanlocs);
X_chan = string(squeeze(chanlocs(1,:,:)));

% Compute absolute channel profiles 
efp_abs_chan = mean(reshape(efp_abs,...
    [n_chans,n_delays*n_bands]),2);

% Compute average channel profiles 
efp_avg_chan = mean(reshape(efp_avg,...
    [n_chans,n_delays*n_bands]),2);

% Plot absolute channel profiles 
figure_name = strcat('Absolute channel profile across all subjects');    
figure('Name',figure_name); topoplot(efp_abs_chan,EEG.chanlocs,...
    'electrodes','labels','whitebk','on', 'gridscale',500);
ax = gca;outerpos = ax.OuterPosition;ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
colorbar; caxis([0 max(efp_abs_chan)]);
img_out = strcat('PROFILE_ABSEFP_Chan');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot average channel profiles 
figure_name = strcat('Average channel profile across all subjects');    
figure('Name',figure_name); topoplot(efp_avg_chan,EEG.chanlocs,...
    'electrodes','labels','whitebk','on', 'gridscale',500);
colorbar; caxis([min(efp_avg_chan) max(efp_avg_chan)]);
img_out = strcat('PROFILE_AVGEFP_Chan');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot std of channel profiles 
efp_chan = reshape(efp,[n_chans,...
    n_delays,n_bands,n_subjects]);
efp_chan = reshape(efp_chan,[n_chans,...
    n_delays*n_bands*n_subjects]);
efp_chan_std = std(efp_chan');
efp_chan_std_abs = std(abs(efp_chan'));

% Plot std of absolute channel profiles 
figure_name = strcat('Std of absolute channel profile across all subjects');    
figure('Name',figure_name); topoplot(efp_chan_std_abs,EEG.chanlocs,...
    'electrodes','labels','whitebk','on', 'gridscale',500);
ax = gca;outerpos = ax.OuterPosition;ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
colorbar; caxis([min(efp_chan_std) max(efp_chan_std)]);
img_out = strcat('PROFILE_ABSEFP_ChanStd');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Plot std of average channel profiles 
figure_name = strcat('Std of channel profile across all subjects');    
figure('Name',figure_name); topoplot(efp_chan_std,EEG.chanlocs,...
    'electrodes','labels','whitebk','on', 'gridscale',500);
ax = gca;outerpos = ax.OuterPosition;ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
colorbar; caxis([min(efp_chan_std) max(efp_chan_std)]);
img_out = strcat('PROFILE_AVGEFP_ChanStd');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Boxplot of absolute channel profiles
efp_chan = reshape(efp,[n_chans,...
    n_delays*n_bands*n_subjects]);
boxplot(abs(efp_chan),'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',X_chan,'FontSize',18)
ylabel('Absolute EFP','FontSize',22);
img_out = strcat('BOXPLOT_ABSEFP_Chan');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

% Boxplot of channel profiles
boxplot(efp_chan,'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',X_chan,'FontSize',18)
ylabel('EFP','FontSize',22);
img_out = strcat('BOXPLOT_AVGEFP_Chan');
saveas(gcf,fullfile(path_out,img_out),'png');
saveas(gcf,fullfile(path_out,img_out),'fig');

end

% ============================================================
% write_group_stats_file()             
% ============================================================

function [] = write_group_stats_file(cfg)

switch cfg.compare
    
    case 'metric'
    in_labels = ["lc4" "lc6" "rmsf" "tp"];
    in_runs = strcat('run_',in_labels,'_',...
        cfg.shift,'_',cfg.cv,'.mat');
    in_paths = strcat(cfg.regress,'/.mat/',in_runs);

    case 'regress'
    in_labels = ["l21_1" "elasticnet"];
    in_runs = strcat('run_',cfg.metric,'_',...
        cfg.shift,'_',cfg.cv,'.mat');
    in_paths = strcat(in_labels,'/.mat/',in_runs);
    
end

id = [in_labels' in_paths'];
file_id = fopen(strcat('RESULTS/',...
    'group_stats/group_stats.txt'),'w');
fprintf(file_id,'%s %s\n', id');
fclose(file_id);

end

