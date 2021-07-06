function [] = report_reliability(data, data_reliab, ...
    metric, R, path_data_out, path_img_out, path_pars) 
% 
%   [] = report_reliability plots and saves reliability
%        results 

% -------------------------------------------------
% Initial settings 
% -------------------------------------------------

import mlreportgen.dom.*;
import mlreportgen.report.*;
    
% Get the parameters for
% the current metric
get_metric_pars;

% Figure settings 
ax = gca;outerpos = ax.OuterPosition; ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Reshape
data = reshape(data, [dim size(data, 2)]);
data_reliab = reshape(data_reliab, dim);

% -------------------------------------------------
% Topographies - report 
% -------------------------------------------------

% Minimum and maximum intensity values 
min_val = min(min(min(data_reliab)));
max_val = max(max(max(data_reliab)));
int = [min_val max_val];

% Try new colormap (higher sensitivity for higher values)
cmap = hot(256);
center_point = 0.5;
scaling_intensity = 5;
x = 1 : length(cmap);
x = x - (center_point - min_val)*length(x)/(max_val - min_val);
x = scaling_intensity*x/max(abs(x));
x = sign(x).*log(abs(x));
x = x - min(x); x = x*511/max(x) + 1;
newmap = interp1(x, cmap, 1:512);


% Default topoplot settings 
topo_settings = {'electrodes', 'on', ...
                'whitebk', 'on', 'gridscale', 100, ...
                'colormap', colormap(newmap)};

% only one band
% only one delay 
% exclude full connectivity
if n_bands == 1 && n_delays == 1 && length(dim) <= 3

    % Create figure
    % current band 
    fig = figure();

    %Increase your figure pixel resolution
    fig.Position(3:4) = fig.Position(3:4)*5;   

    % Specify signal for plotting 
    signal = squeeze(data_reliab(:, 1, 1));
    topoplot(signal, chanlocs, topo_settings{:});
    colorbar; 
    caxis(int);  

     img_out = strcat(upper(metric), '_.png');
     source = fullfile(path_img_out, img_out);
     saveas(fig, source);
     I = FormalImage(source);
     I.ScaleToFit=true; add(R, I); 

% only one band
% many delays
% exclude full connectivity
elseif n_bands == 1 && length(dim) <= 3

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
        signal = squeeze(data_reliab(:, d, 1));

        % Plots for the report 
        subplot(subp_rows, subp_cols, d);
        title(strcat(upper(id_delays(d)), ' SECONDS'));
        topoplot(signal, chanlocs, topo_settings{:})
        colorbar; 
        caxis(int);

    end % finish looping through delays

    img_out = strcat(upper(metric), '_ALLDELAYS_.png');
    source = fullfile(path_img_out, img_out);
    saveas(fig, source); I = FormalImage(source);
    I.ScaleToFit=true; add(R,I);    

% only one delay
% many bands 
% exclude full connectivity
elseif n_delays == 1 && length(dim) <= 3

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
        signal = squeeze(data_reliab(:, 1, b));

        % Plots for the report 
        subplot(subp_rows, subp_cols, b);
        title(upper(id_bands(b)));
        topoplot(signal, chanlocs, topo_settings{:});
        colorbar; 
        caxis(int);  

    end

    img_out = strcat(upper(metric), '_ALLBANDS_.png');
    source = fullfile(path_img_out, img_out);
    saveas(fig, source); I = FormalImage(source);
    I.ScaleToFit=true; add(R, I);

% many delays
% many bands
% exclude full connectivity 
elseif length(dim) <= 3

    % Initialize the band index
    % For models with more then 5,
    % bands, don't plot delta and 
    % theta data values 
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
            signal = squeeze(data_reliab(:, d, b));

            % Plots for the report 
            subplot(subp_rows, subp_cols, d);
            title(strcat(upper(id_delays(d)), ' SECONDS'));
            topoplot(signal, chanlocs, topo_settings{:});
            colorbar;
            caxis(int);  

        end % finish looping through delays

        img_out = strcat(upper(metric), '_', upper(id_bands(b)), '_.png');
        source = fullfile(path_img_out, img_out);
        saveas(fig, source); I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);

    end % finish looping through bands     

end

% -------------------------------------------------
% Topographies - at each band, delay pair 
% -------------------------------------------------    

for d = 1 : n_delays

    % We won't average 
    % tstat topo maps 
    if n_bands > 1
        continue
    end

    if n_delays == 1
        continue
    end
    
    if length(dim) > 3
        continue
    end

    % Specify signal for plotting
    signal = squeeze(data_reliab(:, d));

    figure();
    topoplot(signal, chanlocs, topo_settings{:}); 
    colorbar; 
    caxis(int);    

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', id_delays(d), 'SEC'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');        

end %finish looping through delays 

if n_delays == 1
    
% -------------------------------------------------
% Reliability boxplots - 1 delay
% -------------------------------------------------  

    signal = squeeze(data_reliab);

    % Exclude full connectivity 
    if length(dim) <= 3
        
        box_signal = repmat(signal, 1, 1, size(class_chans, 2));
        for c = 1 : size(class_chans, 2)
           box_signal(~logical(class_chans(:, c)), :, c) = NaN; 
        end
      
    % Full connectivity 
    else
        
        box_signal = repmat(signal, 1, 1, 1, size(class_chans, 2));
        for c = 1 : size(class_chans, 2)
           box_signal(~logical(class_chans(:, c)), :, :, c) = NaN; 
        end
        box_signal = reshape(box_signal,[size(box_signal, 1)*size(box_signal, 2), ...
            size(box_signal, 3), size(class_chans, 2)]);
         
    end
        
    figure();
    box = {squeeze(box_signal(:, :, 1)), squeeze(box_signal(:, :, 2)), ...   
         squeeze(box_signal(:, :, 3)), squeeze(box_signal(:, :, 4))};     

    bo = boxplotGroup(box, 'PrimaryLabels', ...
        id_class_chans, 'SecondaryLabels', cellstr(upper(id_bands)));
    set(bo.axis.Children(1).Children, 'Color', '#0072BD', 'linew', 0.75)
    set(bo.axis.Children(2).Children, 'Color', '#D95319', 'linew', 0.75)
    set(bo.axis.Children(3).Children, 'Color', '#77AC30', 'linew', 0.75)
    set(bo.axis.Children(4).Children, 'Color', '#A2142F', 'linew', 0.75)
    %set(bo.axis.Children(5).Children, 'Color', '#7E2F8E', 'linew', 0.75)
    hold on; ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; hold on;   
    ax.YAxis.FontSize = 18;

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_RELIAB_BOX'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');    

    % 2-way anova   
    % Exclude full connectivity 
    if length(dim) <= 3
        
        g1 = zeros(size(signal)); g2 = g1;
        for c = 1 : size(class_chans, 2)
            g1(logical(class_chans(:, c)), :) = c;
        end
        for b = 1 : n_bands 
            g2(:, b) = b;
        end
        
    % Full connectivity     
    else 
        
        g1 = zeros(size(signal)); g2 = g1;
        for c = 1 : size(class_chans, 2)
            g1(logical(class_chans(:, c)), :, :) = c;
        end
        for b = 1 : n_bands 
            g2(:, :, b) = b;
        end

    end 
    
    [~, ~, signal_stats_1] = anovan(signal(:), {g1(:), g2(:)});
    [~, ~, signal_stats_2] = anovan(signal(:), {g2(:), g1(:)});

    % Compute pairwise results of the multiple comparison
    % test - obtain p-values for the hypothesis test that
    % the corresponding pairwise mean difference is not 0
    % Cols 1 and 2 contain the indices of the two samples
    % being compared; col 3 is the lower confidence interval,
    % col 4 is the estimate, col 5 is the upper confidence 
    % interval, col 6 is the p-value 
    multcomp_1 = multcompare(signal_stats_1); 
    multcomp_2 = multcompare(signal_stats_2); 

    % Save multiple comparison results for current 
    % field in output path 
    save(fullfile(path_data_out, strcat(metric, '_', ...
        'chan_class_reliab_multcomp.mat')), 'multcomp_1');    
    save(fullfile(path_data_out, strcat(metric, '_', ...
        'bands_reliab_multcomp.mat')), 'multcomp_2');

% -------------------------------------------------
% Data boxplots - 1 delay  
% -------------------------------------------------  

    % Exclude full connectivity
    if length(dim) <= 3
        signal = permute(squeeze(data), [2 1 3]);

        box_signal = repmat(signal, 1, 1, 1, size(class_chans, 2));
        for c = 1 : size(class_chans, 2)
           box_signal(:, ~logical(class_chans(:, c)), :, c) = NaN; 
        end

        box_signal = reshape(box_signal, [size(box_signal, 1) ...
            size(box_signal, 2)*size(box_signal, 3) size(box_signal, 4)]);
    
    else
        
        signal = permute(squeeze(data), [3 1 2 4]);

        box_signal = repmat(signal, 1, 1, 1, 1, size(class_chans, 2));
        for c = 1 : size(class_chans, 2)
           box_signal(:, ~logical(class_chans(:, c)), :, :, c) = NaN; 
        end

        box_signal = reshape(box_signal, [size(box_signal, 1) ...
            size(box_signal, 2)*size(box_signal, 3)*size(box_signal, 4) ...
            size(box_signal, 5)]);
        
    end
    
    figure();
     box = {squeeze(box_signal(:, :, 1))', squeeze(box_signal(:, :, 2))', ...   
         squeeze(box_signal(:, :, 3))', squeeze(box_signal(:, :, 4))'}; 

    bo = boxplotGroup(box, 'PrimaryLabels', ...
        id_class_chans, 'SecondaryLabels', cellstr(upper(id_bands)));
    set(bo.axis.Children(1).Children, 'Color', '#0072BD', 'linew', 0.75)
    set(bo.axis.Children(2).Children, 'Color', '#D95319', 'linew', 0.75)
    set(bo.axis.Children(3).Children, 'Color', '#77AC30', 'linew', 0.75)
    set(bo.axis.Children(4).Children, 'Color', '#A2142F', 'linew', 0.75)
    %set(bo.axis.Children(5).Children, 'Color', '#7E2F8E', 'linew', 0.75)
    hold on; ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; hold on;   
    ax.YAxis.FontSize = 18;
    
    % Save figure in specified output path
    img_out = strcat(upper(metric), '_BOX'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');   
    
% -------------------------------------------------
% Dependencies - ICC vs abs(data) - 1 delay
% ------------------------------------------------- 

    % Go through bands 
    for b = 1 : n_bands

        % Exclude full connectivity 
        if length(dim) <= 3
            
            signal = squeeze(data(:, :, b, :));
            signal_reliab = squeeze(data_reliab(:, :, b));

            x = mean(abs(signal), 2); y = signal_reliab;
           
        % Full connectivity     
        else 
            
            signal = squeeze(data(:, :, :, b, :));
            signal_reliab = squeeze(data_reliab(:, :, :, b));
            
            x = mean(abs(signal), 3); y = signal_reliab;
            
            % Reliability can be zero for connectivities 
            % that are zero in all subjects and datasets
            x(isnan(y)) = []; y(isnan(y)) = [];
            
        end
        
        % Fit 
        p = polyfit(x(:), y(:), 1);
        f = polyval(p, x); 

        figure;
        scatter(x(:), y(:), 5, 'filled'); hold on;
        plot(x, f); hold on
        text(max(x(:)), max(y(:)), ...
            strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
        xlabel('ABS(DATA)'); ylabel('ICC');

        % Save figure in specified output path
        img_out = strcat(upper(metric), '_', ...
            upper(id_bands(b)), '_ICC_VS_ABSDATA'); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');  

    end % bands 
    
    % For all bands
    signal = squeeze(data);
    signal_reliab = squeeze(data_reliab);

    % Exclude full connectivity 
    if length(dim) <= 3
        x = mean(abs(signal), 3); y = signal_reliab;
    % Full connectivity 
    else 
        x = mean(abs(signal), 4); y = signal_reliab;
        x(isnan(y)) = []; y(isnan(y)) = [];
    end

    % Fit 
    p = polyfit(x(:), y(:), 1);
    f = polyval(p, x); 

    figure;
    scatter(x(:), y(:), 5, 'filled'); hold on
    plot(x(:), f(:)); hold on
    text(max(x(:)), max(y(:)), ...
        strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
        xlabel('ABS(DATA)'); ylabel('ICC');

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_ICC_VS_ABSDATA'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');     
    
end   

if n_delays == 1; return; end 

% -------------------------------------------------
% Reliability boxplots - many delays 
% -------------------------------------------------  

% Go through frequency bands 
for b = 1 : n_bands
    
    signal = squeeze(data_reliab(:, :, b));
   
    box_signal = repmat(signal, 1, 1, size(class_chans, 2));
    for c = 1 : size(class_chans, 2)
       box_signal(~logical(class_chans(:, c)), :, c) = NaN; 
    end

    figure();
    box = {squeeze(box_signal(:, 1, :)), squeeze(box_signal(:, 2, :)), ...   
         squeeze(box_signal(:, 3, :)), squeeze(box_signal(:, 4, :)), ...   
         squeeze(box_signal(:, 5, :)), squeeze(box_signal(:, 6, :))}; 
    
    bo = boxplotGroup(box, 'PrimaryLabels', ...
        cellstr(upper(id_delays)), 'SecondaryLabels', id_class_chans);
    set(bo.axis.Children(1).Children, 'Color', '#0072BD', 'linew', 1.5)
    set(bo.axis.Children(2).Children, 'Color', '#D95319', 'linew', 1.5)
    set(bo.axis.Children(3).Children, 'Color', '#77AC30', 'linew', 1.5)
    set(bo.axis.Children(4).Children, 'Color', '#A2142F', 'linew', 1.5)
    set(bo.axis.Children(5).Children, 'Color', '#7E2F8E', 'linew', 1.5)
    set(bo.axis.Children(6).Children, 'Color', '#F26FDA', 'linew', 1.5)
    hold on; ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; hold on;   
    ax.YAxis.FontSize = 24;
    ylim(ax,[-0.2 0.6]); set(ax,'ytick', -0.2:0.1:0.6);
    
    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', upper(id_bands(b)), '_RELIAB_BOX'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');    
    
    % 2-way anova 
    g1 = zeros(size(signal)); g2 = g1;
    for c = 1 : size(class_chans, 2)
        g1(logical(class_chans(:, c)), :) = c;
    end
    for d = 1 : n_delays 
        g2(:, d) = d;
    end

    [~, ~, signal_stats_1] = anovan(signal(:), {g1(:), g2(:)});
    [~, ~, signal_stats_2] = anovan(signal(:), {g2(:), g1(:)});
    
    % Compute pairwise results of the multiple comparison
    % test - obtain p-values for the hypothesis test that
    % the corresponding pairwise mean difference is not 0
    % Cols 1 and 2 contain the indices of the two samples
    % being compared; col 3 is the lower confidence interval,
    % col 4 is the estimate, col 5 is the upper confidence 
    % interval, col 6 is the p-value 
    multcomp_1 = multcompare(signal_stats_1); 
    multcomp_2 = multcompare(signal_stats_2); 
    
    % Save multiple comparison results for current 
    % field in output path 
    save(fullfile(path_data_out, strcat(upper(metric), ...
        '_chan_class_reliab_multcomp.mat')), 'multcomp_1');    
    save(fullfile(path_data_out, strcat(upper(metric), ...
        '_delays_reliab_multcomp.mat')), 'multcomp_2');

end

if n_bands > 1
    
    % 3-way anova for differences between bands 
    g1 = zeros(size(data_reliab)); g2 = g1; g3 = g1;

    for b = 1 : n_bands
        g1(:, :, b) = b;
    end
    for c = 1 : size(class_chans, 2)
        g2(logical(class_chans(:, c)), :, :) = c;
    end
    for d = 1 : n_delays 
        g3(:, d, :) = d;
    end

    [~, ~, signal_stats] = anovan(data_reliab(:), {g1(:), g2(:), g3(:)});
    multcomp = multcompare(signal_stats); 
     save(fullfile(path_data_out, strcat(upper(metric), ...
         '_band_reliab_multcomp.mat')), 'multcomp');    

end

% -------------------------------------------------
% Data boxplots - many delays 
% -------------------------------------------------  

% Go through frequency bands 
for b = 1 : n_bands
    
    signal = permute(squeeze(data(:, :, b, :)), [2 1 3]);
    
    box_signal = repmat(signal, 1, 1, 1, size(class_chans, 2));
    for c = 1 : size(class_chans, 2)
       box_signal(:, ~logical(class_chans(:, c)), :, c) = NaN; 
    end
    
    box_signal = reshape(box_signal, [size(box_signal, 1) ...
        size(box_signal, 2)*size(box_signal, 3) size(box_signal, 4)]);
    
    figure();
     box = {squeeze(box_signal(1, :, :)), squeeze(box_signal(2, :, :)), ...   
         squeeze(box_signal(3, :, :)), squeeze(box_signal(4, :, :)), ...   
         squeeze(box_signal(5, :, :)), squeeze(box_signal(6, :, :))}; 
    
    bo = boxplotGroup(box, 'PrimaryLabels', ...
        cellstr(upper(id_delays)), 'SecondaryLabels', id_class_chans);
    set(bo.axis.Children(1).Children, 'Color', '#0072BD', 'linew', 0.75)
    set(bo.axis.Children(2).Children, 'Color', '#D95319', 'linew', 0.75)
    set(bo.axis.Children(3).Children, 'Color', '#77AC30', 'linew', 0.75)
    set(bo.axis.Children(4).Children, 'Color', '#A2142F', 'linew', 0.75)
    set(bo.axis.Children(5).Children, 'Color', '#7E2F8E', 'linew', 0.75)
    set(bo.axis.Children(6).Children, 'Color', '#F26FDA', 'linew', 0.75)
    hold on; ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; hold on;   
    ax.YAxis.FontSize = 18;
    
    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', upper(id_bands(b)), '_BOX'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');    

end % bands  

% -------------------------------------------------
% Dependencies - ICC vs abs(data) - many delays 
% ------------------------------------------------- 

% Go through bands 
for b = 1 : n_bands
    
    signal = squeeze(data(:, :, b, :));
    signal_reliab = squeeze(data_reliab(:, :, b));

    x = mean(abs(signal), 3); y = signal_reliab;

    % Fit 
    p = polyfit(x(:), y(:), 1);
    f = polyval(p, x); 
        
    figure;
    scatter(x(:), y(:), 5, 'filled'); hold on
    plot(x, f); hold on
	text(max(x(:)), max(y(:)), ...
        strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
        xlabel('ABS(DATA)'); ylabel('ICC');

    % Save figure in specified output path
    img_out = strcat(upper(metric), '_', ...
        upper(id_bands(b)), '_ICC_VS_ABSDATA'); 
    saveas(gcf, fullfile(path_img_out, img_out), 'png');  

end % bands

% For all bands
signal = data;
signal_reliab = data_reliab;

x = mean(abs(signal), 4); y = signal_reliab;

% Fit 
p = polyfit(x(:), y(:), 1);
f = polyval(p, x); 

figure;
scatter(x(:), y(:), 5, 'filled'); hold on
plot(x(:), f(:)); hold on
text(max(x(:)), max(y(:)), ...
    strcat('CORR: ', " ", num2str(corr(x(:), y(:)))));
    xlabel('ABS(DATA)'); ylabel('ICC');

% Save figure in specified output path
img_out = strcat(upper(metric), '_ICC_VS_ABSDATA'); 
saveas(gcf, fullfile(path_img_out, img_out), 'png');  

% -------------------------------------------------
% Topographies - at each band, delay pair 
% -------------------------------------------------    

 set(groot, 'defaultFigureUnits','centimeters');
 set(groot, 'defaultFigurePosition',[0 0 18 30]);
     
% Only for cases where delas are defined 
if n_bands == 1; return; end
if n_delays == 1; return; end

for d = 1 : n_delays
    
    for b = 1 : n_bands  
        
        % Specify signal for plotting 
        signal = squeeze(data_reliab(:, d, b));
        
        figure();
        topoplot(signal, chanlocs, topo_settings{:}); 
        colorbar; 
        caxis(int);    

        % Save figure in specified output path
        img_out = strcat(upper(metric), '_', upper(id_bands(b)), '_', ...
            upper(id_delays(d)), 'SEC'); 
        saveas(gcf, fullfile(path_img_out, img_out), 'png');
        
    end % finish looping through bands
    
end % finish looping through delays

end