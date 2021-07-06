% Performs reliability correlation analysis, from the results
% of the subject-level correlation analysis, in which the Pearson's
% correlation between the EEG features and the corresponding BOLD
% signal was estimated 
% Plots the results, saves them in the IMAGES directory and  
% generates an automatic report 

% Import report APIs 
import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
        
% Create a title for the report
if flag.report ~= 0
    
    my_title = 'CORRELATION RELIABILITY ANALYSIS';
    H1 = get_report_heading(1, my_title);
    add(R, H1)  
    
end

n_subjects = length(subjects);
n_sessions = length(sessions); 
n_metrics = length(metrics);

% ------------------------------------------------------------
% Retreive subjects correlation stats
% ------------------------------------------------------------    

subj_stats = cell(n_metrics, 1);
subj_power = cell(n_metrics, 1); 
subj_power_reliab = cell(n_metrics, 1);
                                 
% Retreive subject-level
% correlation stats 
for m = 1 : n_metrics
   
    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);
       
    % Pre-allocate correlation matrix
    rho_all = zeros(n_features, n_subjects, n_sessions);
    
    % Pre-allocate power matrix 
    power_all = zeros(n_features/n_delays, n_subjects, n_sessions); 
    
    % Pre-allocate power reliability matrix
    power_reliab_all = zeros(n_features/n_delays, 1);
    
    for se = 1 : n_sessions
        
        for s = 1 : n_subjects

            subject = subjects(s);

            % Load subject correlation data 
            load(fullfile(path_data_in(s, se), ...
                strcat(metric, '_', data_in)));
            rho_all(:, s, se) = stats.rho; 
            
            % Load subject power or connectivity data
             if contains(metric, power_metrics) 
                load(fullfile(path_data_in_power(s, se), ...
                    strcat(metric, '_', data_in_power)));
                power_all(:, s, se) = power(:);
             elseif contains(metric, connectivity_metrics)
                load(fullfile(path_data_in_connectivity(s, se), ...
                    strcat(metric, '_', data_in_connectivity)));
                power_all(:, s, se) = net(:); 
             end
            
        end % finish looping through subjects
        
    end % finish looping through sessions 
    
    % Load subject power or connectivity reliability data
    if contains(metric, power_metrics) 
        load(fullfile(path_data_in_power_reliab, ...
            strcat(metric, '_', data_in_power_reliab)));
        power_reliab_all(:) = reliab(:); 
    elseif contains(metric, connectivity_metrics)
        load(fullfile(path_data_in_connectivity_reliab, ...
            strcat(metric, '_', data_in_connectivity_reliab)));
        power_reliab_all(:) = reliab(:); 
    end    
      
    subj_stats{m, 1} = reshape(rho_all, [n_features n_subjects*n_sessions]);
    subj_power{m, 1} = reshape(power_all, [n_features/n_delays n_subjects*n_sessions]);
    subj_power_reliab{m, 1} = power_reliab_all;
    
end % finish looping through metrics 

% Go through metrics 
for m = 1 : n_metrics

    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);  
        
    % Broadcast the current pipeline stage 
    disp(strcat('Performing correlation reliability', ...
        ' analysis for metric'," ", metric, ' ...'));
        
    % Add metric heading 
    % to the report 
    if flag.report ~=0
        
        my_title = upper(metric);
        H2 = get_report_heading(2, my_title);
        add(R, H2);
        
    end   
    
    % ------------------------------------------------------
    % Reliability tests tests 
    % ------------------------------------------------------    
    
    % Compute measures of intra- and inter-subject variability,
    % as well as the interclass correlation coefficient for each 
    % EEG feature 
    reliab = reliability_test(cell2mat(subj_stats(m, ...
        1)), n_subjects, n_sessions, reliab_metric, ...
        fullfile(path_data_out, strcat(metric, ...
        '_', data_out)));
    
    half_corr = reliability_split_half_test(cell2mat(subj_stats(m, 1)), ...
        metric, n_subjects, n_sessions, ...
        n_iter_reliab, path_pars, path_data_out, path_img_out);     
    
    % ------------------------------------------------------
    % Plot reliability results   
    % ------------------------------------------------------    
    
    if flag.report ~= 0
        
        % Add report title 
        my_title = upper(reliab_metric);
        H3 = get_report_heading(3, my_title);
        add(R, H3);
        
        report_reliability(cell2mat(subj_stats(m, 1)), ...
            reliab, metric, R, path_data_out, path_img_out, path_pars); 
        
    end 
    
    % ------------------------------------------------------
    % Plot dependency results   
    % ------------------------------------------------------ 

      power = mean(mean(reshape(cell2mat(subj_power(m, 1)), [n_chans ...
          n_bands n_subjects n_sessions]), 4), 3);
      power_reliab = reshape(cell2mat(subj_power_reliab(m, 1)), ...
          [n_chans n_bands]);
      reliab = squeeze(mean(reshape(reliab, [n_chans n_delays n_bands]), 2));
      report_dependencies(power, power_reliab, reliab, metric, ...
          path_img_out, path_pars); 
    
end % finish looping through metrics


function [] = report_dependencies(data_x, data_x_reliab, ...
    data_y_reliab, metric, path_img_out, path_pars)

% Figure settings 
ax = gca;outerpos = ax.OuterPosition; ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

get_metric_pars

% -------------------------------------------------
% Dependencies - Corr. Reliab. vs. Avg. Power 
% ------------------------------------------------- 

% Go through bands 
for b = 1 : n_bands
    
    x = squeeze(data_x(:, b));
    y = squeeze(data_y_reliab(:, b));

   % Fit 
    p = polyfit(x(:), y(:), 1);
    f = polyval(p, x); 
        
    figure;
    scatter(x(:), y(:), 5, 'filled'); hold on 
    plot(x, f); hold on
	text(max(x(:)), max(y(:)), ...
        strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
        xlabel('Power/Net. Measure'); ylabel('ICC');

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', ...
        upper(id_bands(b)), '_CORR_RELIAB_VS_POWER'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');  

end % bands

% All bands 
x = data_x;
y = data_y_reliab;

% Fit 
p = polyfit(x(:), y(:), 1);
f = polyval(p, x); 

figure;
scatter(x(:), y(:), 5, 'filled'); hold on 
plot(x, f); hold on
text(max(x(:)), max(y(:)), ...
    strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
    xlabel('Power/Net. Measure'); ylabel('ICC');

% Save figure in specified output path
img_out = strcat(upper(metric), '_CORR_RELIAB_VS_POWER'); 
saveas(gcf, fullfile(path_img_out, img_out), 'png'); 

% -------------------------------------------------
% Dependencies - Corr. Reliab. vs. Power Reliab. 
% ------------------------------------------------- 

% Go through bands 
for b = 1 : n_bands
    
    x = squeeze(data_x_reliab(:, b));
    y = squeeze(data_y_reliab(:, b));

   % Fit 
    p = polyfit(x(:), y(:), 1);
    f = polyval(p, x); 
            
    figure;
    scatter(x(:), y(:), 5, 'filled'); hold on
    plot(x, f); hold on
	text(max(x(:)), max(y(:)), ...
        strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
        xlabel('Power/Net. Measure Reliability'); ylabel('ICC');    

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', ...
        upper(id_bands(b)), '_CORR_RELIAB_VS_POWER_RELIAB'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');  

end % bands

% All bands 
x = data_x_reliab;
y = data_y_reliab;

% Fit 
p = polyfit(x(:), y(:), 1);
f = polyval(p, x); 

figure;
scatter(x(:), y(:), 5, 'filled'); hold on
plot(x, f); hold on
text(max(x(:)), max(y(:)), ...
    strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
    xlabel('Power/Net. Measure Reliability'); ylabel('ICC');    

% Save figure in specified output path
img_out = strcat(upper(metric), '_CORR_RELIAB_VS_POWER_RELIAB'); 
saveas(gcf, fullfile(path_img_out, img_out), 'png');  

end