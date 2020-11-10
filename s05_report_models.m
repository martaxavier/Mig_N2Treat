% Create the report heading 
my_title = strcat(model,' MODELS');
H1 = get_report_heading(1,my_title);
add(R,H1) 

n_subjects = length(subjects);
n_metrics = length(metrics);
fs = fs_analysis;

% ------------------------------------------------------------
% Go through metrics 
% ------------------------------------------------------------

for m = 1 : n_metrics
    
    metric = metrics(m);
    
    % Get parameters
    % for current metric
    get_metric_pars
    
    % Specify directory where images are to be saved 
    path_img_metric_out = strcat(path_img_out,'\',metric);

    % Create directory where results are to be saved 
    if ~exist(path_img_metric_out, 'dir')
        mkdir(path_img_metric_out); 
    end
    
    % Load matrix containing the optimal 
    % delay for the BOLD deconvolution 
    % of each subject
    if strcmp(bold_shift,'deconv')
        load(fullfile(path_pars_in, ...
            'deconv_delay.mat'));
    end
    
    % Add metric heading to the report 
    my_title = upper(metric);
    H2 = get_report_heading(2,my_title);
    add(R,H2);
    
    % --------------------------------------------------------
    % Go through subjects 
    % --------------------------------------------------------

    for s = 1 : n_subjects
        
        subject = subjects(s);
        
        % Broadcast the current pipeline stage 
        disp(strcat('Creating report of model ', ...
    ' results for subject'," ",subject,', metric', ...
    " ",metric, ' ...'));
        
        % Define input EEG and BOLD data, according
        % to current metric
        eeg_in = strcat(eeg_metric,'_', ...
            'eeg_feature','_',eeg_shift,'.txt');
        bold_in = strcat('bold_processed', ...
        '_',bold_shift,'.txt');
        if contains(eeg_in,'_.'); eeg_in = ...
                replace(eeg_in,'_.','.'); end
        if contains(bold_in,'_.'); bold_in = ...
                replace(bold_in,'_.','.'); end
    
        % Load BOLD and EEG data
        eeg = dlmread(fullfile(path_eeg_in(s),eeg_in));
        bold = dlmread(fullfile(path_bold_in(s),bold_in));

        % Load model results for current subject,
        % model pair 
        model_in = strcat(metric,'_','model.mat');
        load(fullfile(path_model_in(s,m),model_in));
        
        % Save estimated deconvolution 
        % delay for current subject 
        % (when applicable)
        deconv_delay = [];
        if strcmp(bold_shift,'deconv')
            deconv_delay = ...
                table2array(deconv_delay(s,2));
        end

        % Add subject heading to the report 
        my_title = subjects(s);
        H3 = get_report_heading(3,my_title);
        add(R,H3);

        % Plot and/or report model results (EFP)
        plot_efp(bold,model,metric, ...
            deconv_delay,fs,R,flag.report, ...
            path_img_metric_out);

    end % finish looping through metrics 

end % finish looping through subjects


%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

function [] = plot_efp(bold,model,reg_model,metric, ...
    deconv_delay,fs,R,report,path_img_out)

% ------------------------------------------------------------
% Read inputs 
% ------------------------------------------------------------  

% Import report APIs 
import mlreportgen.dom.*;
import mlreportgen.report.*;
  
% Get parameters for 
% current metric 
get_metric_pars;

% Figure settings 
ax = gca; outerpos = ax.OuterPosition; 
ti = ax.TightInset.*3.5; 
left = outerpos(1) + ti(1); 
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3); 
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Write default settings for topoplot
topo_settings = {'electrodes','labels', ...
    'whitebk','on','gridscale',300};

% ------------------------------------------------------------
% Assign model variables 
% ------------------------------------------------------------ 

% Assign model variables 
% that are averaged through folds 
bic_model = model.bic;
nmse_model = model.nmse;
corr_model = model.corr;
lambda_model = model.lambda;
rho_model = model.rho; 
df_model = model.df;
efp_model = model.efp;
yhat_model = model.yhat;

% ------------------------------------------------------------
% Plot model variables - performance and final model
% ------------------------------------------------------------

% Plot model results and save figure in output path  

if strcmp (reg_model,'l21_1'); id_rho = 'Rho';
else; id_rho = 'Alpha'; end

rho_msg = strcat(id_rho,' final model:'," ", num2str(rho_model));
lambda_msg = strcat('Lambda final model:'," ", num2str(lambda_model));
df_msg = strcat('DOF final model:'," ", num2str(df_model));
bic_msg = strcat('Average BIC:'," ", num2str(bic_model)); 
nmse_msg = strcat('Average NMSE:'," ", num2str(nmse_model));
corr_msg = strcat('Average correlation:'," ", num2str(corr_model));
message = {rho_msg, lambda_msg, df_msg, bic_msg, nmse_msg, corr_msg};

figure('Name','Model variables'); title('Model variables')
x = [0.1 0.1 0.1 0.1 0.1 0.1]; y = [0.8 0.7 0.6 0.5 0.4 0.3];
text(x,y,message, 'FontSize', 14); set(gca,'xtick',[],'ytick',[]);
img_out = 'MODEL.png'; saveas(gcf,fullfile(path_img_out,img_out));

% ------------------------------------------------------------
% Report model variables - performance and final model
% ------------------------------------------------------------

% Report the model performance 
my_title = 'MODEL PERFORMANCE';
H4 = get_report_heading(4,my_title);
add(R,H4);

P = Paragraph();
append(P, Text(sprintf(' ')));
P.WhiteSpace = 'preserve';
add(R,P);

% T = Text(strcat('Model performance was', ...
%     ' assessed by a 15-fold nondependent', ...
%     ' cross-validation procedure. The model', ...
%     'performance parameters are the following:'));
% T.FontSize = "10pt";
% T.Color = '#404040';
% T.FontFamilyName = 'Arial';
% add(R,T);
% 
% P = Paragraph();
% append(P, Text(sprintf(' ')));
% P.WhiteSpace = 'preserve';
% add(R,P);

txt = [bic_msg,nmse_msg,corr_msg];
for t = 1 : length(txt)
    T = Text(strcat('-'," ",txt(t)));
    T.FontSize = "12pt";
    T.Color = '#404040';
    T.FontFamilyName = 'Arial';
    add(R,T);
end

P = Paragraph();
append(P, Text(sprintf(' ')));
P.WhiteSpace = 'preserve';
add(R,P);

% Report the final model 
my_title = 'FINAL PERFORMANCE';
H4 = get_report_heading(4,my_title);
add(R,H4);


P = Paragraph();
append(P, Text(sprintf(' ')));
P.WhiteSpace = 'preserve';
add(R,P);

% T = Text(strcat('The final model was estimated', ...
%     ' in the entire dataset, using the final', ...
%     '  regularization parameters. Bellow are the', ...
%     ' final regularization parameters, estimated through', ...
%     ' a nested cross-validation procedure, and the degrees', ...
%     ' of freedom (DOF) of the final model estimated:'));
% T.FontSize = "10pt";
% T.Color = '#404040';
% T.FontFamilyName = 'Arial';
% add(R,T);
% 
% P = Paragraph();
% append(P, Text(sprintf(' ')));
% P.WhiteSpace = 'preserve';
% add(R,P);

txt = [rho_msg,lambda_msg, df_msg];
for t = 1 : length(txt)
    T = Text(strcat('-'," ",txt(t)));
    T.FontSize = "12pt";
    T.Color = '#404040';
    T.FontFamilyName = 'Arial';
    add(R,T);
end

% ------------------------------------------------------------
% Plot BOLD signal and BOLD estimate 
% ------------------------------------------------------------

% Only if image
% generation is on 
if report == 1 
    
    time = 0 : 1/fs : (length(bold)-1)*(1/fs);

    % Plot predicted BOLD signal on top of original BOLD signal
    fig = figure('Name', strcat('Predicted BOLD signal', ...
        ' and original BOLD signal'));
    fig.Position(3:4) = fig.Position(3:4)*5;
    nmse_msg = upper(strcat('nmse: ', num2str(nmse_model)));
    text(time(end)-50,-2,nmse_msg, 'FontSize', 22, ...
        'FontName','Arial'); hold on;
    plot(time, bold,'LineWidth',2.5,'Color','#1ca9df'); hold on; 
    plot(time, yhat_model,'LineWidth',3,'Color','#4acb84');
    legend('BOLD signal', 'BOLD estimate','FontSize',22, ...
        'FontName','Arial');
    xlabel('Time (s)','FontSize',22,'FontName','Arial');
    ylabel('Amplitude','FontSize',22,'FontName','Arial');
    grid on 
    img_out = 'BOLD.png';
    saveas(gcf,fullfile(path_img_out,img_out));

end

% ------------------------------------------------------------
% Report BOLD signal and BOLD estimate
% ------------------------------------------------------------

% Create title for current report section  
my_title = 'BOLD ESTIMATE';
H4 = get_report_heading(4,my_title);
add(R,H4);

I = Image(fullfile(path_img_out,img_out));
I.Style={ScaleToFit(true),HAlign('center')};
add(R,I);
        
% ------------------------------------------------------------
% Plot EEG coefficient estimates 
% ------------------------------------------------------------

% Create title for current report section  
my_title = 'ALL COEFFICIENTS';
H4 = get_report_heading(4,my_title);
add(R,H4);

% Remove first coefficient (mean)
efp_model = efp_model(2:end);

% Reshape EEG coefficients matrix (EFP)
efp_model = reshape(efp_model,efp_dim);

% Reshape EEG coefficients matrix 
% for better visualization for plotting 
efp_plot = permute(efp_model,[2,3,1]);
efp_plot = reshape(efp_plot,[n_delays*...
    n_bands,n_chans]);

xtickname = cellstr(id_chans(chans));
xangle = 0; xname = 'CHANNELS';

% Plot EEG coefficient estimates matrix (EFP)
fig = figure('Name','Coefficient estimates matrix'); 
fig.Position(3:4) = fig.Position(3:4)*5;
imagesc(efp_plot); colorbar;
xlabel(xname,'FontSize',20); 
ylabel('BANDS x DELAYS','FontSize',20);
set(gca,'ytick',[]); 
xticks(1:n_chans); 
xticklabels(xtickname);
xtickangle(xangle); 

% Save and add to report 
img_out = 'EFP.png';
saveas(gcf,fullfile(path_img_out,img_out));
I = Image(fullfile(path_img_out,img_out));
I.Style={ScaleToFit(true),HAlign('center')};
add(R,I);

% Plot absolute EEG coefficient estimates matrix (EFP_ABS)
figure('Name','Coefficient estimates matrix (abs)'); 
imagesc(abs(efp_plot)); colorbar;
xlabel(xname,'FontSize',20); 
ylabel('BANDS x DELAYS','FontSize',20);
set(gca,'ytick',[]);
xticks(1:n_chans);
xticklabels(xtickname);

% Save and add to report 
xtickangle(xangle); img_out = 'EFP_abs.png';
saveas(gcf,fullfile(path_img_out,img_out));
    
% ---------------------------------------------
% Channel profiles - at each band, each delay
% ---------------------------------------------

% Define the same maximum and minimum 
% intensities for all band-delay plots 
max_signal_abs = max(max(max(max(efp_model))), ...
    abs(min(min(min(efp_model)))));

% Define default title for following figures
fig_title = strcat('Topographic map of', ...
    'model coefficients for');

% Topoplot of coefficient estimates,
% for each band, each delay, each channel
for b = 1 : n_bands
    
    % Create title for current band 
    my_title = strcat(upper(id_bands...
        (b)),' COEFFICIENTS');
    H4 = get_report_heading(4,my_title);
    add(R,H4);

    % Plot average coefficient estimate
    % values for current band 
    signal = squeeze(efp_model(:,:,b));
    
    figure('Name',strcat(fig_title, ...
        " ",id_bands(b),' power')); 
    fig.Position(3:4) = fig.Position(3:4)*5;

    for d = 1 : n_delays 

        subplot(2,3,d); title(strcat(id_delays(d), ...
            " ", 'SECONDS'),'FontSize',14);
        topoplot(signal(:,d),chanlocs,topo_settings{:});
        colorbar; caxis([-max_signal_abs max_signal_abs]); 

    end

    % Save figure and add to report 
    img_out = strcat('TOPO_AVGEFP_',id_bands(b),'.png'); 
    saveas(gcf,fullfile(path_img_out,img_out));
    I = Image(fullfile(path_img_out,img_out));
    I.Style={ScaleToFit(true),HAlign('center')};
    add(R,I);
        
    % Plot absolute coefficient estimate
    % values for current band 
    signal = squeeze(abs(efp_model(:,:,b)));

    figure('Name',strcat(fig_title," ", ...
        id_bands(b),' power (abs)')); 

    for d = 1 : n_delays 

        subplot(2,3,d); title(strcat(id_delays(d), ...
            " ", 'SECONDS'),'FontSize',14);
        topoplot(signal(:,d),chanlocs,topo_settings{:});
        colorbar; caxis([0 max_signal_abs]);

    end

    % Save figure 
    img_out = strcat('TOPO_ABSEFP_', ...
        id_bands(b),'_abs.png'); 
    saveas(gcf,fullfile(path_img_out,img_out));

end

% Leave if image 
% generation is off 
if report == 2
    return
end

% ------------------------------------------------------------
% Channel profiles of coefficient estimates 
% ------------------------------------------------------------

% ---------------------------------------------
% Channel profiles - all delays, bands 
% ---------------------------------------------  

% Topographic maps of averaged (through delays and frequencies)
% average coefficient estimates  
fig_title = strcat('Topographic map of averaged',...
   ' (through delays and bands) coefficient estimates'); 
figure('Name',fig_title); signal = sum(sum(efp_model,3),2); 
topoplot(signal,chanlocs,topo_settings{:}); 
colorbar; img_out = strcat('TOPO_AVGEFP.png');
saveas(gcf,fullfile(path_img_out,img_out));

% Topographic maps of averaged (through delays and frequencies)
% average coefficient estimates (absolute value)
fig_title = strcat('Topographic map (abs values) of averaged', ...
   ' (through delays and bands) coefficient estimates'); 
figure('Name',fig_title); signal = abs(sum(sum(efp_model,3),2));
topoplot(signal,chanlocs,topo_settings{:}); 
caxis([0 max(signal)]); colorbar; 
img_out = strcat('TOPO_AVGEFP_abs.png');
saveas(gcf,fullfile(path_img_out,img_out));

% Topographic maps of averaged (through delays and frequencies)
% absolute coefficient estimates  
fig_title = strcat('Topographic map of averaged (through', ...
   'delays and bands) absolute coefficient estimates'); 
figure('Name',fig_title); signal = sum(sum(abs(efp_model),3),2);
topoplot(signal,chanlocs,topo_settings{:});
colorbar; caxis([0 max(signal)]); 
img_out = strcat('TOPO_ABSEFP.png');
saveas(gcf,fullfile(path_img_out,img_out));

% ---------------------------------------------
% Channel profiles - at each delay, all bands 
% ---------------------------------------------

  % Topoplot of coefficient estimates, for each
% delay, averaged through all bands
signal = squeeze(sum(efp_model,3));
max_signal_abs = max(max(max(signal)),abs(min(min(signal))));

% Plot average coefficient estimate values for each delay  
fig_title = strcat('Topographic map of model coefficients',...
    ', averaged through all bands, for each delay') ; 
figure('Name',fig_title);

for d = 1 : n_delays  

    subplot(2,3,d); 
    title(strcat(id_delays(d)," ", 'SECONDS'),'FontSize',14);
    topoplot(signal(:,d),chanlocs,topo_settings{:});
    colorbar; caxis([-max_signal_abs max_signal_abs]);             

end   
img_out = strcat('TOPO_AVGEFP_Delays','.png');
saveas(gcf,fullfile(path_img_out,img_out));               

% Topoplot of absolute coefficient estimates,
% for each delay, averaged through all bands
signal = squeeze(sum(abs(efp_model),3));
max_signal_abs = max(max(signal));

% Plot absolute coefficient estimate values for each delay  
fig_title = strcat('Topographic map of absolute model',...
    ' coefficients, averaged through all bands, for each delay') ; 
figure('Name',fig_title);

for d = 1 : n_delays

    subplot(2,3,d); 
    title(strcat(id_delays(d)," ", 'SECONDS'),'FontSize',14);
    topoplot(signal(:,d),topo_settings{:});
    colorbar; caxis([0 max_signal_abs]);   

end
img_out = strcat('TOPO_ABSEFP_Delays','.png');
saveas(gcf,fullfile(path_img_out,img_out));


% ---------------------------------------------
% Channel profiles - at each band, all delays
% ---------------------------------------------

if n_bands > 1

    % Topoplot of coefficient estimates,
    % for each band, averaged through all delays
    signal = squeeze(sum(efp_model,2));
    max_signal_abs = max(max(max(signal)),abs(min(min(signal))));

    % Plot average coefficient estimate values for each band 
    fig_title = strcat('Topographic map of model coefficients',...
        ', averaged through all delays, for each frequency band') ; 
    figure('Name',fig_title);

    for b = 1 : n_bands  

        subplot(3,2,b); 
        title(strcat(upper(id_bands(b)),' BAND'),'FontSize',14);
        topoplot(signal(:,b),chanlocs,topo_settings{:});
        colorbar; caxis([-max_signal_abs max_signal_abs]);           

    end   
    img_out = strcat('TOPO_AVGEFP_Bands','.png');
    saveas(gcf,fullfile(path_img_out,img_out));               

    % Topoplot of coefficient estimates,
    % for each band, averaged through all delays
    signal = squeeze(sum(abs(efp_model),2));
    max_signal_abs = max(max(signal));

    % Plot absolute coefficient estimate values for each band 
    fig_title = strcat('Topographic map of absolute model', ...
        ' coefficients, averaged through all delays, for each band') ; 
    figure('Name',fig_title);

    for b = 1 : n_bands 

        subplot(3,2,b); 
        title(strcat(upper(id_bands(b)), ' BAND'),'FontSize',14);
        topoplot(abs(signal(:,b)),topo_settings{:});
        colorbar; caxis([0 max_signal_abs]);   

    end
    img_out = strcat('TOPO_ABSEFP_Bands','.png');
    saveas(gcf,fullfile(path_img_out,img_out));         

end     

% ------------------------------------------------------------
% Frequency profiles of coefficient estimates 
% ------------------------------------------------------------

% Delay labels for bar plots 
X_del = categorical...
(cellstr(upper(id_delays)));
X_del = reordercats...
(X_del,cellstr(upper(id_delays)));

if n_bands > 1 

    % Frequency labels for bar plots
    X_freq = categorical(cellstr(upper(id_bands)));
    X_freq = reordercats(X_freq,cellstr(upper(id_bands)));

    % ---------------------------------------------
    % Frequency profiles - all delays and chans
    % ---------------------------------------------

    % Compute abs and avg frequency profiles 
    efp_freq = reshape(efp_model,...
        [n_chans*n_delays,n_bands]);
    efp_avg_freq = mean(efp_freq);
    efp_abs_freq = mean(abs(efp_freq));

    % Plot frequency profiles of the values of the EFP, 
    % averaged through channels and delays 
    figure_name = 'Average frequency profile'; figure('Name',figure_name);
    bar(X_freq,efp_avg_freq,'FaceColor','#77AC30','BarWidth',0.6);
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
    img_out = strcat('PROFILE_AVGEFP_Bands.png');
    saveas(gcf,fullfile(path_img_out,img_out));

    % Plot frequency profiles of the absolute values of the EFP, 
    % averaged through channels and delays 
    figure_name = 'Absolute frequency profile'; figure('Name',figure_name);
    bar(X_freq,efp_abs_freq,'FaceColor','#77AC30','BarWidth',0.6);
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
    img_out = strcat('PROFILE_ABSEFP_Bands.png');
    saveas(gcf,fullfile(path_img_out,img_out));

    % Boxplot of frequency profiles of the values
    % of the EFP, averaged through channels and delays
    boxplot(efp_freq,'Colors',[0 0.4470 0.7410])
    set(gca,'XTickLabel',cellstr(X_freq),'FontSize',18)
    ylabel('EFP','FontSize',22);
    img_out = strcat('BOXPLOT_AVGEFP_Bands.png');
    saveas(gcf,fullfile(path_img_out,img_out));

    % Boxplot of frequency profiles of the absolute values
    % of the EFP, averaged through channels and delays
    boxplot(abs(efp_freq),'Colors',[0 0.4470 0.7410])
    set(gca,'XTickLabel',cellstr(X_freq),'FontSize',18)
    ylabel('Absolute EFP','FontSize',22);
    img_out = strcat('BOXPLOT_ABSEFP_Bands.png');
    saveas(gcf,fullfile(path_img_out,img_out));

    % ---------------------------------------------
    % Frequency profiles - at each delay, all chans
    % ---------------------------------------------

    % Compute absolute frequency profiles at each delay
    efp_abs_freqbydel = ...
        squeeze(mean(reshape(abs(efp_model),...
        efp_dim)));

    % Compute average frequency profiles at each delay
    efp_avg_freqbydel = ...
        squeeze(mean(reshape(efp_model,...
        efp_dim)));

    % Plot frequency profiles of the values of the EFP, 
    % at each delay, averaged through channels 
    figure_name = 'Average frequency profile for each delay';    
    figure('Name',figure_name); 
    b = bar(X_del,efp_avg_freqbydel,'FaceColor','flat','BarWidth',0.9);
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
    for k = 1:size(efp_avg_freqbydel,2);b(k).CData = k; end
    set(b,{'DisplayName'},cellstr(X_freq)'); legend('FontSize',16);
    img_out = strcat('PROFILE_AVGEFP_BandbyDelay.png');
    saveas(gcf,fullfile(path_img_out,img_out));

    % Plot frequency profiles of the absolute values of the EFP, 
    % at each delay, averaged through channels 
    figure_name = 'Absolute frequency profile for each delay';    
    figure('Name',figure_name); 
    b = bar(X_del,efp_abs_freqbydel,'FaceColor','flat','BarWidth',0.9);
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
    for k = 1:size(efp_avg_freqbydel,2);b(k).CData = k; end
    set(b,{'DisplayName'},cellstr(X_freq)'); legend('FontSize',16);
    img_out = strcat('PROFILE_ABSEFP_BandbyDelay.png');
    saveas(gcf,fullfile(path_img_out,img_out));

end

% ------------------------------------------------------------
% Delay profiles of coefficient estimates 
% ------------------------------------------------------------

% ---------------------------------------------
% Delay profiles - all bands and chans
% ---------------------------------------------

% Compute abs and avg delay profiles 
efp_del = reshape(permute(efp_model, ...
[1 3 2]),[n_chans*n_bands,n_delays]);
efp_abs_del = mean(abs(efp_del));
efp_avg_del = mean(efp_del);

% Plot delay profiles of the values of the EFP, 
% averaged through channels and bands 
figure_name = 'Average delay profile'; figure('Name',figure_name);
bar(X_del,efp_avg_del,'FaceColor','flat','BarWidth',0.6);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
img_out = strcat('PROFILE_AVGEFP_Delay.png');
saveas(gcf,fullfile(path_img_out,img_out));

% Plot delay profiles of the absolute values of the EFP, 
% averaged through channels and bands 
figure_name = 'Absolute delay profile'; figure('Name',figure_name);
bar(X_del,efp_abs_del,'FaceColor','flat','BarWidth',0.6);
ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
img_out = strcat('PROFILE_ABSEFP_Delay.png');
saveas(gcf,fullfile(path_img_out,img_out));

% Boxplot of delay profiles of the average values
% of the EFP, averaged through channels and bands
boxplot(efp_del,'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',cellstr(X_del),'FontSize',18)
ylabel('EFP','FontSize',22);
img_out = strcat('BOXPLOT_AVGEFP_Delay.png');
saveas(gcf,fullfile(path_img_out,img_out));

% Boxplot of delay profiles of the absolute values
% of the EFP, averaged through channels and bands
boxplot(abs(efp_del),'Colors',[0 0.4470 0.7410])
set(gca,'XTickLabel',cellstr(X_del),'FontSize',18)
ylabel('Absolute EFP','FontSize',22);
img_out = strcat('BOXPLOT_ABSEFP_Delay.png');
saveas(gcf,fullfile(path_img_out,img_out));

if n_bands > 1 

    % ---------------------------------------------
    % Delay profiles - at each band, all chans 
    % ---------------------------------------------

    % Compute average frequency profiles at each delay
    efp_avg_delbyfreq = squeeze(mean...
        (reshape(efp_model,efp_dim)))';

    % Compute absolute frequency profiles at each delay
    efp_abs_delbyfreq = squeeze(mean...
        (reshape(abs(efp_model),efp_dim)))';

    % Plot average frequency profiles at each delay
    figure_name = strcat('Average delay profile across all', ...
        'subjects, for each frequency'); figure('Name',figure_name); 
    b = bar(X_freq,efp_avg_delbyfreq,'FaceColor','flat','BarWidth',0.9);
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
    for k = 1:size(efp_avg_delbyfreq,2);b(k).CData = k; end
    set(b,{'DisplayName'},cellstr(X_del)'); legend('FontSize',16);
    img_out = strcat('PROFILE_AVGEFP_DelaybyBand.png');
    saveas(gcf,fullfile(path_img_out,img_out));   

    % Plot absolute frequency profiles at each delay
    figure_name = strcat('Absolute delay profile across all', ...
        'subjects, for each frequency band'); figure('Name',figure_name); 
    b = bar(X_freq,efp_abs_delbyfreq,'FaceColor','flat','BarWidth',0.9);
    ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; set(gca,'FontSize',16)
    for k = 1:size(efp_abs_delbyfreq,2);b(k).CData = k; end
    set(b,{'DisplayName'},cellstr(X_del)'); legend('FontSize',16);
    img_out = strcat('PROFILE_ABSEFP_DelaybyBand.png');
    saveas(gcf,fullfile(path_img_out,img_out));

end 

end % finish print_efp

