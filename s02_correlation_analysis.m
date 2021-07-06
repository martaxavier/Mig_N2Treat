import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
    
% ------------------------------------------------------------
% Subject-level analysis
% ------------------------------------------------------------ 

% Create a title for the report
if flag.report ~=0

    my_title = strcat('CORRELATION ANALYSIS -', ...
        " ", session);
    H1 = get_report_heading(1, my_title);
    add(R,H1)
    
end

n_subjects = length(subjects);
n_metrics = length(metrics);

% ------------------------------------------------------------
% Go through metrics
% ------------------------------------------------------------

subj_stats = cell(n_metrics, 2); % first column rho
                                % second column pval
                                
for m = 1 : n_metrics

    metric = metrics(m);
    get_metric_pars

    % Define input EEG and BOLD data, according
    % to current metric
    eeg_in = strcat(eeg_metric, '_', ...
        'eeg_feature', '_', eeg_shift, '.txt');
    bold_in = strcat('bold_preproc', ...
    '_', bold_shift, '.txt');
    if contains(eeg_in, '_.'); eeg_in = ...
            replace(eeg_in, '_.', '.'); end
    if contains(bold_in, '_.'); bold_in = ...
            replace(bold_in, '_.', '.'); end
    
    % Load matrix containing the optimal 
    % delay for the BOLD deconvolution 
    % of each subject
    if strcmp(bold_shift, 'deconv')
        load(fullfile(path_pars, ...
            'deconv_delay.mat'));
    end
    
    n_features = prod(dim);

    % Add the metric title to the report 
    if flag.report ~= 0
        
        my_title = upper(metric);
        H2 = get_report_heading(2, my_title);
        add(R, H2)

    end

    %--------------------------------------------------------    
    % Go through subjects 
    %--------------------------------------------------------
    
    for s = 1 : n_subjects

        subject = subjects(s);
        
        disp(strcat('Performing correlation analysis for', ...
            " ", subject, ', metric', " ", metric, ' ...'));
        
        % Create output directories for 
        % current subject, if non existent
        if ~exist(path_data_out(s, se), 'dir')
            mkdir(path_data_out(s, se));
        end
        if ~exist(path_img_out(s, se), 'dir')
            mkdir(path_img_out(s, se));
        end      

        % Add subject title to the report 
        if flag.report ~= 0 

            my_title = subject;
            H3 = get_report_heading(3, my_title);
            add(R, H3)
        
        end

        % Load input EEG and BOLD data 
        eeg = dlmread(fullfile(path_eeg_in(s, se), eeg_in));
        bold = dlmread(fullfile(path_bold_in(s, se), bold_in));
        
        % Save estimated deconvolution 
        % delay for current subject 
        % (when applicable)
        deconv_delay = [];
        if strcmp(bold_shift, 'deconv')
            deconv_delay = ...
                table2array(deconv_delay(s, 2));
        end
        
        % Perform correlation analysis for the 
        % current subject 
        stats = perform_corr_analysis(eeg, bold);
        
        % Save results of correlation analysis for 
        % for the current subject and metric 
        corr_out = strcat(metric, '_', data_out);
        save(fullfile(path_data_out(s, se), corr_out), 'stats');     

        % Generate and save plots of the result
        % of correlation analysis
         if flag.report ~= 0
             
             plot_corr_stats(stats.rho, metric, R, ...
                 deconv_delay, flag.report, ...
                 path_img_out(s, se), path_pars, subject);
             
         end

    end % finish looping through subjects 

end % finish looping through metrics  


%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% perform_corr_analysis(eeg,bold)             
% ============================================================

function [stats] = perform_corr_analysis(EEG,BOLD)

%   [stats] = perform_corr_analysis performs correlation 
%             analysis between each of the EEG features 
%             specified and the BOLD signal 
%
%   Inputs:
%
%     EEG              The model's X 
%     BOLD             The model's Y
%
%   Outputs:
%
%   stats              Correlation statistics for current subject 

% -------------------------------------------------
% Process input parameters 
% -------------------------------------------------

% If Y is a row vector, 
% convert to a column vector
if size(BOLD, 1) == 1
    BOLD = BOLD';
end

% -------------------------------------------------
% Compute pairwise correlation  
% -------------------------------------------------

% Compute Pearson's correlation between 
% each pair of columns in EEG and BOLD 
% rho is the pairwise Pearson's correlation
% coefficient; pval is a vector of p-values
% for testing the hypothesis of no correlation
% against the alternative hypothesis of a
% nonzero correlation

[rho, pval] = corr(EEG, BOLD);

% Prepare output
stats.rho = rho;
stats.pval = pval; 

end

% ============================================================
% plot_corr_stats(rho,metric,R,deconv_delay,report,path_img_out)          
% ============================================================

function [] = plot_corr_stats(rho, metric, R, deconv_delay, ...
    report, path_img_out, path_pars, subject)
%
%   [] = plot_corr_stats(rho,metric,R,deconv_delay,report,
%   path_img_out) plots and saves subject-level correlation
%   profiles, in the form of topographic maps or bar graphs 
%
%   Inputs:
%   
%           rho            the correlation values 
%           metric         the current EEG-BOLD metric 
%           R              the report object 
%           report         the report flag - 2, 1 or 0 
%           path_img_out   the path where plots are to be saved 
%

import mlreportgen.dom.*;
import mlreportgen.report.*;
    
% Get metric pars 
get_metric_pars;

% Reshape rho into feature
% dimension space
rho = reshape(rho, dim);

% Figure settings 
ax = gca;outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Write default settings for topoplot
topo_settings = {'electrodes', 'labels', ...
    'whitebk', 'on', 'gridscale', 100};

% Write default title for figures
fig_title = strcat('Correlation values between', ...
' EEG features and the BOLD signal');

% Write default topoplot
% image tag
topo_out = 'TOPOCORR';

if report == 1
    
    % -------------------------------------------------
    % Channel profiles (at all delays, bands)  
    % -------------------------------------------------

    % Specify signals for plotting 
    signal = squeeze(mean(reshape(rho, ...
        [n_chans,n_delays*n_bands]), 2));
    signal_abs = squeeze(mean(reshape(abs(rho), ...
        [n_chans,n_delays*n_bands]), 2));

    % Minimum and maximum intensity values 
    min_val = min(signal);
    min_int = min_val - 0.05*abs(min_val);
    max_val = max(signal);
    max_int = max_val + 0.05*abs(max_val);

    % Topoplot of the correlation between EEG features and the
    % BOLD signal, averaged at each channel location 
    figure('Name', strcat(fig_title, 'averaged at each channel')); 
    topoplot(signal, chanlocs, topo_settings{:}); 
    colorbar; caxis([min_int max_int]);    

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', topo_out);
    saveas(gcf,fullfile(path_img_out, img_out), 'png');

    % Maximum intensity values
    max_abs_val = max(signal_abs);
    max_abs_int = max_abs_val + 0.05*abs(max_abs_val);

    % Topoplot of the absolute value of the correlation
    % between EEG features and the BOLD signal, averaged
    figure('Name',strcat('Absolute'," ", fig_title, ...
        ' averaged at each channel')); 
    topoplot(signal_abs, chanlocs, topo_settings{:}); 
    colorbar; caxis([0 max_abs_int]);    

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', topo_out, '_ABS'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');

    % -------------------------------------------------
    % Channel profiles (at all delays for each band)  
    % -------------------------------------------------

    % Minimum and maximum intensity values 
    % Choose limits for all delays 
    min_val = min(min(squeeze(mean(rho,2))));
    min_int = min_val - 0.1*abs(min_val);
    max_val = max(max(squeeze(mean(rho,2))));
    max_int = max_val + 0.1*abs(max_val);

    max_abs_val = max(max(squeeze(mean(abs(rho), 2))));
    max_abs_int = max_abs_val + 0.1*abs(max_abs_val);

    for b = 1 : n_bands

        if n_bands == 1
            continue
        end

        % Specify signals for plotting 
        signal = squeeze(mean(rho(:, :, b), 2));
        signal_abs = squeeze(mean(abs(rho(:, :, b)), 2));

        % Topoplot of the correlation between EEG features 
        % and the BOLD signal, averaged at each channel, 
        % at each band 
        figure('Name',strcat(fig_title, ...
            ' averaged at each channel, at each band'));
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([min_int max_int]);    

        % Save figure in specified output path
        img_out = strcat(id_bands(b), bold_shift, '_', topo_out); 
        saveas(gcf,fullfile(path_img_out, img_out), 'png');

        % Topoplot of the absolute value of the correlation
        % between EEG features and the BOLD signal, averaged
        % at each channel location, at each band  
        figure('Name', strcat('Absolute', " ", fig_title, ...
            ' averaged at each channel, at each band'));
        topoplot(signal_abs, chanlocs, topo_settings{:}); 
        colorbar; caxis([0 max_abs_int]);    

        % Save figure in specified output path
        img_out = strcat(id_bands(b), bold_shift, ...
            '_', topo_out, '_ABS'); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

    end % finish looping through bands 

    % -------------------------------------------------
    % Channel profiles (at all bands for each delay)  
    % -------------------------------------------------

    % Minimum and maximum intensity values 
    % Choose limits for all delays 
    min_val = min(min(squeeze(mean(rho, 3))));
    min_int = min_val - 0.1*abs(min_val);
    max_val = max(max(squeeze(mean(rho, 3))));
    max_int = max_val + 0.1*abs(max_val);

    max_abs_val = max(max(squeeze(mean(abs(rho), 3))));
    max_abs_int = max_abs_val + 0.1*abs(max_abs_val);

    for d = 1 : n_delays

        if n_delays == 1
            continue
        end

        % Specify signals for plotting 
        signal = squeeze(mean(squeeze(rho(:, d, :)), 2));
        signal_abs = squeeze(mean(squeeze(abs(rho(:, d, :))), 2));

        % Topoplot of the correlation between EEG features 
        % and the BOLD signal, averaged at each channel,
        % at each delay 
        figure('Name',strcat(fig_title, ...
            ' averaged at each channel, at each delay')); 
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([min_int max_int]); 

        % Save figure in specified output path
        img_out = strcat(upper(metric), '_', ...
            id_delays(d), 'SEC_', topo_out); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

        if n_bands == 1
            continue
        end

        % Topoplot of the absolute value of the correlation
        % between EEG features and the BOLD signal, averaged
        % at each channel location, at each delay  
        figure('Name',strcat('Absolute', " ", fig_title, ...
            ' averaged at each channel, at each delay'));
        topoplot(signal_abs, chanlocs, topo_settings{:}); 
        colorbar; caxis([0 max_abs_int]);       

        % Save figure in specified output path
        img_out = strcat(upper(metric), '_', ...
            id_delays(d), 'SEC_', topo_out, '_ABS'); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

    end %finish looping through delays 

end

% -------------------------------------------------
% Channel profiles - report   
% -------------------------------------------------

% Minimum and maximum intensity values 
min_val = min(min(min(rho)));
min_int = min_val - 0.1*abs(min_val);
max_val = max(max(max(rho)));
max_int = max_val + 0.1*abs(max_val);  

% only one band
% only one delay 
if n_bands == 1 && n_delays == 1
    
    img_out = strcat(upper(metric),'_', topo_out, '.png');
    source = fullfile(path_img_out, img_out);
    I = FormalImage(source);
    I.ScaleToFit=true; 
    
    if contains(metric,'deconv')
        caps = strcat('Time-to-Peak of estimated',...
            ' HRF:', " ", num2str(deconv_delay), ' seconds.');
        I.Caption = caps;
    end
    add(R,I); 
     
    
% only one band
% many delays
elseif n_bands == 1
    
    % For the report 
    subp_rows = 2;
    subp_cols = n_delays/subp_rows;
    
    % Create title for current band 
    H4 = get_report_heading(4, id_bands);
    add(R, H4);

    % Create figure
    % current band 
    fig = figure();

    %Increase your figure pixel resolution
    fig.Position(3:4) = fig.Position(3:4)*5;

    for d = 1 : n_delays

        % Specify signal for plotting 
        signal = squeeze(rho(:, d, 1));

        % Plots for the report 
        subplot(subp_rows, subp_cols, d);
        title(strcat(id_delays(d), ' SECONDS'));
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([min_int max_int]);  

    end % finish looping through delays

    img_out = strcat(upper(metric), '_', ...
        'alldelays_', topo_out, '.png');
    source = fullfile(path_img_out, img_out);
    saveas(fig, source); I = FormalImage(source);
    I.ScaleToFit=true; add(R, I);    
    
% only one delay
% many bands 
elseif n_delays == 1
    
    % For the report 
    subp_rows = 2;
    subp_cols = n_bands/subp_rows;
    
    % Create figure
    % current band 
    fig = figure();

    %Increase your figure pixel resolution
    fig.Position(3:4) = fig.Position(3:4)*5;

    for b = 1 : n_bands
        
        % Specify signal for plotting 
        signal = squeeze(rho(:, 1, b));

        % Plots for the report 
        subplot(subp_rows, subp_cols, b);
        title(upper(id_bands(b)));
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([min_int max_int]); 
        
    end
    
    img_out = strcat(upper(metric), '_', ...
        'allbands_', topo_out, '.png');
    source = fullfile(path_img_out, img_out);
    saveas(fig, source); I = FormalImage(source);
        
    if contains(metric,'deconv')
        caps = strcat('Time-to-Peak of estimated', ...
            ' HRF:', " ", num2str(deconv_delay), ' seconds.');
        I.Caption = caps;
    end
    
    I.ScaleToFit=true; add(R, I);
   
% many delays
% many bands 
else
              
    % Initialize the band index
    % For models with more then 4,
    % bands, don't plot delta and 
    % theta correlations 
    if n_bands > 5; b1 = 3;
    else; b1 = 1; end
    
    % For the report 
    subp_rows = 2;
    subp_cols = n_delays/subp_rows;

    for b = b1 : n_bands

        % Create title for current band 
        H4 = get_report_heading(4, upper(id_bands(b)));
        add(R, H4);

        % Create figure
        % current band 
        fig = figure();

        %Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;

        for d = 1 : n_delays

            % Specify signal for plotting 
            signal = squeeze(rho(:, d, b));

            % Plots for the report 
            subplot(subp_rows, subp_cols, d);
            title(strcat(id_delays(d), ' SECONDS'));
            topoplot(signal, chanlocs, topo_settings{:}); 
            colorbar; caxis([min_int max_int]);  

        end % finish looping through delays

        img_out = strcat(upper(metric), '_', ...
            id_bands(b), '_', topo_out, '.png');
        source = fullfile(path_img_out, img_out);
        saveas(fig, source); I = FormalImage(source);
        I.ScaleToFit=true; add(R, I);

    end % finish looping through bands     
    
end

% Leave if image 
% generation is off 
if report == 1
    return
end
    
% -------------------------------------------------
% Channel profiles (at each band-delay pair)  
% -------------------------------------------------

% Only for cases where bands are defined 
if n_bands == 1; return; end
if n_delays == 1; return; end

% Minimum and maximum intensity values 
min_val = min(min(min(rho)));
min_int = min_val - 0.1*abs(min_val);
max_val = max(max(max(rho)));
max_int = max_val + 0.1*abs(max_val);  

for d = 1 : n_delays
    
    for b = 1 : n_bands  
        
        % Specify signal for plotting 
        signal = squeeze(rho(:, d, b));
        
        % Topoplot of the correlation between EEG features
        % and the BOLD signal, averaged at each channel,
        % at each delay-band pair  
        figure('Name', strcat(fig_title, ...
            ' at each channel, at each delay and band')); 
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([min_int max_int]);    

        % Save figure in specified output path
        img_out = strcat(id_bands(b), '_', ...
            id_delays(d), 'SEC_', topo_out); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');
    
    end % finish looping through bands
    
end % finish looping through delays

end
