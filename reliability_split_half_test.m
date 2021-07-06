function [half_corr] = reliability_split_half_test(data, metric, ...
    n_subjects, n_sessions, n_iter, path_pars, path_data_out, path_img_out)

%-----------------------------------------------------------
% Check population consistency through split half test
%-----------------------------------------------------------

% PLOT BOXPLOTS 

get_metric_pars;
data = reshape(data, [dim, n_subjects, n_sessions]);

% Excluding full connectivity
if length(dim) <= 3
    data = permute(data, [3 1 2 4 5]);
    data = reshape(data, [size(data, 1) ...
        size(data, 2) size(data, 3) n_subjects n_sessions]);
% For full connectivity
elseif length(dim) > 3
    data = permute(data, [4 1 2 3 5 6]);
    data = reshape(data, [size(data, 1) ...
        size(data, 2)*size(data, 3) ...
        size(data, 4) n_subjects n_sessions]);
end

half_corr = zeros(4, 4, n_iter);

for i = 1 : n_iter
    
    split_half = crossvalind('KFold', n_subjects, 2);
    
    data_first_half = squeeze(mean(mean(data(:, :, :, split_half == 1, :), 5), 4));
    data_second_half = squeeze(mean(mean(data(:, :, :, split_half == 2, :), 5), 4));
    
    data_first_half = reshape(data_first_half, [size(data, 1) size(data, 2)*size(data, 3)]);
    data_second_half = reshape(data_second_half, [size(data, 1) size(data, 2)*size(data, 3)]);
        
    half_corr(:, :, i) = corr(data_first_half', data_second_half'); 
    
end

% Figure settings 
ax = gca;outerpos = ax.OuterPosition; ti = ax.TightInset; 
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Boxplot of results 
figure();
box = {abs(squeeze(half_corr(1, :, :)))', abs(squeeze(half_corr(2, :, :)))', ...   
     abs(squeeze(half_corr(3, :, :)))', abs(squeeze(half_corr(4, :, :)))'}; 

bo = boxplotGroup(box, 'PrimaryLabels', ...
    cellstr(upper(id_bands)), 'SecondaryLabels', cellstr(upper(id_bands)));
set(bo.axis.Children(1).Children, 'Color', '#0072BD', 'linew', 0.75)
set(bo.axis.Children(2).Children, 'Color', '#D95319', 'linew', 0.75)
set(bo.axis.Children(3).Children, 'Color', '#77AC30', 'linew', 0.75)
set(bo.axis.Children(4).Children, 'Color', '#A2142F', 'linew', 0.75)
hold on; ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; hold on;   
ax.YAxis.FontSize = 18;
    
% Save figure in specified output path
img_out = strcat(upper(metric), '_SPLIT_HALF_CORR_BOX'); 
saveas(gcf, fullfile(path_img_out, img_out), 'png');    
    
half_corr = mean(half_corr, 3);

data_out = strcat(metric, '_', 'split_half_corr.mat');
save(fullfile(path_data_out, data_out), 'half_corr');
