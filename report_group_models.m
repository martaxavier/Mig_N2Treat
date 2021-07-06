function [] = report_group_models(group_efp, metric, ...
    R, prob, pval, path_img_out, path_pars)

%
%   [] = report_group_models(group_models, metric, ...
%        R, prob, report, path_img_out) plots and saves the 
%        group-level tstat profiles (thresholded and
%        unthresholded) in the form of topographic maps
%        or bar graphs 
%
%   Inputs:
%   
%           group_models   the group correlation/model values  
%           metric         the current EEG-BOLD metric
%           R              the report object 
%           prob           the probability of topographic 
%                          consistency
%           pval           the p-value of topographic consistency
%           path_img_out   path where plots are to be saved 
%

import mlreportgen.dom.*;
import mlreportgen.report.*;
    
% Get the parameters for
% the current metric
get_metric_pars;

% -------------------------------------------------
% Compute measures  
% -------------------------------------------------
group_efp = group_efp';
efp_avg = mean(group_efp);
efp_std = std(group_efp);

% Figure settings 
ax = gca;outerpos = ax.OuterPosition;ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Write default settings for topoplot
topo_settings = {'electrodes', 'labels', ...
                'conv', 'on', 'whitebk', ...
                'on', 'gridscale', 100};
            
% -------------------------------------------------
% Channel profiles - report   
% -------------------------------------------------

% Minimum and maximum intensity values 
min_val = min(min(min(efp_avg)));
min_int = min_val;
max_val = max(max(max(efp_avg)));
max_int = max_val;
int = max(abs(min_int), abs(max_int));
 
% Write report title for each t
my_title = ["Average Coefficient Estimates", ...
          "Standard Deviation of Coefficient Estimates"];

% Write topoplot tags for each t 
topo_out = ["EFP_AVG", "EFP_STD"];
      
for t = 1 : 2
    
    if t == 1
        efp = squeeze(efp_avg);
    elseif t == 2
        efp = squeeze(efp_std);
    end
    efp = reshape(efp, dim);
    
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
        signal = squeeze(efp(:, 1, 1));
        topoplot(signal, chanlocs, topo_settings{:});
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
            signal = squeeze(efp(:, d, 1));

            % Plots for the report 
            subplot(subp_rows, subp_cols, d);
            title(strcat(id_delays(d), ' SECONDS'));
            topoplot(signal, chanlocs, topo_settings{:});
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
            signal = squeeze(efp(:, 1, b));

            % Plots for the report 
            subplot(subp_rows, subp_cols, b);
            title(upper(id_bands(b)));
            topoplot(signal, chanlocs, topo_settings{:});
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
                signal = squeeze(efp(:, d, b));

                % Plots for the report 
                subplot(subp_rows, subp_cols, d);
                title(strcat(id_delays(d), ' SECONDS'));
                topoplot(signal, chanlocs, topo_settings{:});
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

end % finish looping through types of maps           