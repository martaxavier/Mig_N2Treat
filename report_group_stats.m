function [] = report_group_stats(group_stats, thresh, metric, ...
    R, prob, pval, report, path_img_out)
%
%   [] = report_group_stats(group_stats,thresh,metric, ...
%        R,prob,report,path_img_out) plots and saves the 
%        group-level tstat profiles (thresholded and
%        unthresholded) in the form of topographic maps
%        or bar graphs 
%
%   Inputs:
%   
%           group_stats    the group correlation/model values  
%           thresh         the tstat threshold 
%           metric         the current EEG-BOLD metric
%           R              the report object 
%           prob           the probability of topographic 
%                          consistency
%           pval           the p-value of topographic consistency
%           report         the report flag - 2, 1 or 0
%           path_img_out   path where plots are to be saved 
%

import mlreportgen.dom.*;
import mlreportgen.report.*;
    
% Get the parameters for
% the current metric
get_metric_pars;

% Reshape relevant statistics into feature space 
tstat = reshape(group_stats(:,1), dim);
decision = reshape(group_stats(:,2), dim);

% Threshold tstat
tstat_thresh = tstat;
tstat_thresh(abs(tstat_thresh) < thresh)=0;

% Figure settings 
ax = gca;outerpos = ax.OuterPosition;ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Write default settings for topoplot
topo_settings = {'electrodes', 'labels', ...
                'conv', 'on', 'whitebk', ...
                'on', 'gridscale', 300};
topo_thresh_settings = cat(2,topo_settings, ...
                       {'maplimits', 'absmax'});

% Write default title for figures
fig_title = strcat('tstat map of group', ...
    'statistics between each EEG', ...
    'feature and the BOLD signal');

% Write default topoplot tag
 topo_out = 'TOPOSTAT';

 if report == 1
     
    % -------------------------------------------------
    % Channel profiles (at all bands for each delay)  
    % -------------------------------------------------

    % NOTE - when plotting pvalues in totoplot, set 
    %        'conv' option to 'on' to minimize 
    %        interpolation effects 

    % Minimum and maximum intensity values 
    % Choose limits for all delays 
    min_val = min(min(squeeze(mean(tstat, 3))));
    min_int = min_val;
    max_val = max(max(squeeze(mean(tstat, 3))));
    max_int = max_val;
    int = max(abs(min_int),abs(max_int));

    min_val_thresh = min(min(squeeze(mean(tstat_thresh, 3))));
    min_int_thresh = min_val_thresh;
    max_val_thresh = max(max(squeeze(mean(tstat_thresh, 3))));
    max_int_thresh = max_val_thresh;
    int_thresh = max(abs(min_int_thresh), abs(max_int_thresh));

    for d = 1 : n_delays

        % We won't average 
        % tstat topo maps 
        if n_bands > 1
            continue
        end

        if n_delays == 1
            continue
        end

        % Specify signal for plotting
        signal = squeeze(tstat(:, d));
        signal_thresh = squeeze(tstat_thresh(:, d));
        
        % Topoplot of the unthresholded group one-sample
        % tstat for the EEG-BOLD correlation/model values
        % at each channel, delay and all bands  
        figure('Name', fig_title);
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([-int int]);    

        % Save figure in specified output path
        img_out = strcat(metric, '_', id_delays(d), ...
            'SEC_', topo_out, '_UNTHRESH'); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

        % Topoplot of the thresholded group one-sample
        % tstat for the EEG-BOLD correlation/model values
        % at each channel, delay and all bands  
        figure('Name', strcat('Thresholded', " ", fig_title));
        topoplot(signal_thresh, chanlocs, topo_thresh_settings{:});
        colorbar; caxis([-int_thresh int_thresh]);    

        % Save figure in specified output path
        img_out = strcat(metric, '_', id_delays(d), ...
            'SEC_', topo_out, 'THRESH', num2str(thresh)); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

        % Topoplot of the group one-sample tstat, masked 
        % with the decision map at p = 0.05 significance
        % level, for the EEG-BOLD correlation/model values
        % at each channel, delay and all bands  
        figure('Name', strcat('Significant (p<0.05)', ...
            " ", fig_title));
        topoplot(signal, chanlocs, 'pmask', ...
            decision(:,d,1), topo_settings{:}); 
        colorbar; caxis([-int_thresh int_thresh]);

        % Save figure in specified output path
        img_out = strcat(metric, '_', id_delays(d), ...
            'SEC_', topo_out, '_SIGNIFICANT', num2str(0.05)); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

    end %finish looping through delays 

 end
 
% -------------------------------------------------
% Channel profiles - report   
% -------------------------------------------------

% Minimum and maximum intensity values 
min_val = min(min(min(tstat)));
min_int = min_val;
max_val = max(max(max(tstat)));
max_int = max_val;
int = max(abs(min_int), abs(max_int));
 
% Write report title for each t
my_title = ["Tstat Map, Unthresholded", ...
          "Tstat Map, Thresholded (p<0.05)"];

% Write topoplot tags for each t 
topo_out = ["TOPOSTAT_UNTHRESH",...
    "TOPOSTAT_SIGNIFICANT"];
      
for t = 1 : 2
    
    H3 = get_report_heading(3, my_title(t));
    add(R,H3);
    
    % only one band
    % only one delay 
    if n_bands == 1 && n_delays == 1

        % Create figure
        % current band 
        fig = figure();

        %Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;   
        
        % Specify signal for plotting 
        signal = squeeze(tstat(:, 1, 1));
        
        if t == 1
            topoplot(signal, chanlocs, topo_settings{:});
        elseif t == 2
            topoplot(signal, chanlocs, 'pmask', ...
                decision(:, 1, 1), topo_settings{:}); 
        end
        colorbar; caxis([-int int]);  
        
         img_out = strcat(metric,'_',...
             topo_out(t), '.png');
         source = fullfile(path_img_out, img_out);
         saveas(fig, source);
         I = FormalImage(source);
         I.ScaleToFit=true; add(R, I); 

    % only one band
    % many delays
    elseif n_bands == 1

        % For the report 
        subp_rows = 2;
        subp_cols = n_delays/subp_rows;

        % Create title for current band 
        H4 = get_report_heading(4, upper(id_bands));
        add(R, H4);

        % Create figure
        % current band 
        fig = figure();

        %Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;

        for d = 1 : n_delays

            % Specify signal for plotting 
            signal = squeeze(tstat(:, d, 1));

            % Plots for the report 
            subplot(subp_rows, subp_cols, d);
            title(strcat(id_delays(d), ' SECONDS'));

            if t == 1
                topoplot(signal, chanlocs, topo_settings{:});
            elseif t == 2
                topoplot(signal, chanlocs, 'pmask', ...
                    decision(:,d,1), topo_settings{:}); 
            end
            colorbar; caxis([-int int]);            

        end % finish looping through delays

        % Generate caption string
        my_cap = "";
        for d = 1 : n_delays
            my_cap = strjoin(cat(2, my_cap, strcat(" ", ...
                id_delays(d), ' seconds -', " ", ...
                num2str(prob(d,1)),'%; (p-value =', ...
                " ", num2str(pval(d,1)), ')')));
        end
        my_cap = strcat('Probability of',...
            ' significance:', my_cap);
        
        img_out = strcat(metric, '_', ...
            'alldelays_', topo_out(t), '.png');
        source = fullfile(path_img_out, img_out);
        saveas(fig, source); I = FormalImage(source);
        I.ScaleToFit=true; 
        I.Caption = my_cap; 
        add(R,I);    

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
            signal = squeeze(tstat(:, 1, b));

            % Plots for the report 
            subplot(subp_rows, subp_cols, b);
            title(upper(id_bands(b)));
            
            if t == 1
                topoplot(signal, chanlocs, topo_settings{:});
            elseif t == 2
                topoplot(signal, chanlocs, 'pmask', ...
                    squeeze(decision(:,1,b)), topo_settings{:}); 
            end
            colorbar; caxis([-int int]);  

        end

        % Generate caption strings
        my_cap ="";
        for b = 1 : n_bands
            my_cap=strjoin(cat(2,my_cap, strcat(" ", ...
                id_bands(b), ' band -'," ", ...
                num2str(prob(1,b)), '% (p-value = ;', ...
                " ", num2str(pval(1,b)), ')')));
        end
        my_cap = strcat('Probability of',...
            ' significance:', my_cap);
            
        img_out = strcat(metric, '_',...
            'allbands_', topo_out(t), '.png');
        source = fullfile(path_img_out, img_out);
        saveas(fig, source); I = FormalImage(source);
        I.ScaleToFit=true; 
        I.Caption = my_cap; 
        add(R, I);

    % many delays
    % many bands 
    else

        % Initialize the band index
        % For models with more then 4,
        % bands, don't plot delta and 
        % theta correlations 
        if n_bands > 4; b1 = 3;
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
                signal = squeeze(tstat(:, d, b));

                % Plots for the report 
                subplot(subp_rows, subp_cols, d);
                title(strcat(id_delays(d), ' SECONDS'));

                if t == 1
                    topoplot(signal, chanlocs, topo_settings{:});
                elseif t == 2
                    topoplot(signal, chanlocs, 'pmask', ...
                        decision(:, d, b), topo_settings{:}); 
                end
                colorbar; caxis([-int int]);  

            end % finish looping through delays

            % Generate caption strings
            my_cap ="";
            for d = 1 : n_delays
                my_cap=strjoin(cat(2, my_cap, strcat(" ", ...
                    id_delays(d), ' seconds -', " ", ...
                    num2str(prob(d,b)), '%; (p-value =', ...
                    " ", num2str(pval(d,b)), ')')));
            end
            my_cap = strcat('Probability of',...
                ' significance:', my_cap);

            img_out = strcat(metric, '_', ...
                id_bands(b), '_', topo_out(t), '.png');
            source = fullfile(path_img_out, img_out);
            saveas(fig, source); I = FormalImage(source);
            I.ScaleToFit=true; 
            I.Caption = my_cap; 
            add(R,I);

        end % finish looping through bands     

    end

end % finish looping through tstat maps 

% Leave if flag for 
% generating images 
% is off 
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
min_val = min(min(min(tstat)));
min_int = min_val;
max_val = max(max(max(tstat)));
max_int = max_val;
int = max(abs(min_int),abs(max_int));

min_val_thresh = min(min(min(tstat_thresh)));
min_int_thresh = min_val - 0.1*abs(min_val_thresh);
max_val_thresh = max(max(max(tstat_thresh)));
max_int_thresh = max_val + 0.1*abs(max_val_thresh);
int_thresh = max(abs(min_int_thresh), abs(max_int_thresh));
 
% Write default topoplot tag
topo_out = 'TOPOSTAT';
 
for d = 1 : n_delays
    
    for b = 1 : n_bands  
        
        % Specify signal for plotting 
        signal = squeeze(tstat(:, d, b));
        signal_thresh = squeeze(tstat_thresh(:, d, b));
        
        % Topoplot of the unthresholded group one-sample
        % tstat for the EEG-BOLD correlation/model values 
        % at each channel, delay and band 
        figure('Name', fig_title);
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; caxis([-int int]);    

        % Save figure in specified output path
        img_out = strcat(id_bands(b), '_', ...
            id_delays(d), 'SEC_', topo_out, '_UNTHRESH'); 
        saveas(gcf,fullfile(path_img_out, img_out), 'png');

        % Topoplot of the thresholded group one-sample 
        % tstat for the EEG-BOLD correlation/model values
        % at each channel, delay and band
        figure('Name', strcat('Thresholded', " ", fig_title));
        topoplot(signal_thresh, chanlocs, topo_thresh_settings{:});
        colorbar; caxis([-int_thresh int_thresh]);    

        % Save figure in specified output path
        img_out = strcat(id_bands(b), '_', id_delays(d),...
            'SEC_', topo_out, '_THRESH', num2str(thresh)); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');

        % Topoplot of the group one-sample tstat, masked
        % with the decision map at p = 0.05 significance
        % level, for the EEG-BOLD correlation/model values
        % at each channel, delay and band 
        figure('Name', strcat('Significant (p<0.05)', " ", fig_title));
        topoplot(signal, chanlocs, 'pmask', ...
            decision(:, d, b), topo_settings{:});
        colorbar; caxis([-int_thresh int_thresh]);

        % Save figure in specified output path
        img_out = strcat(id_bands(b), '_', ...
            id_delays(d), 'SEC_', topo_out, '_SIGNIFICANT'); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');
    
    end % finish looping through bands
    
end % finish looping through delays

end
