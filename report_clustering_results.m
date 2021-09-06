function [] = report_clustering_results(ses_coef, hc_distance, path_img_out)

if ~exist(path_img_out, 'dir')
    mkdir(path_img_out);
end

n_sessions = size(ses_coef, 1); 

% Create a hierarchical cluster tree, Z
% Uses euclidean distance as the distance between leaves
% Uses the centroid of each cluster to compute distance 
% between clusters 
Z = linkage(ses_coef, hc_distance(1), hc_distance(2));

% Find the optimal number of clusters by finding the elbow 
% (or knee) method - find the knee of the objective 
% function, which in this case is the correlation function 
[n_clusters, ~] = knee_pt(Z(:, 3));

% Define clusters from the hierarchical cluster tree Z
T = cluster(Z, 'MaxClust', n_clusters);

% Find the largest cluster
[occurences, clusters]= hist(T, unique(T));
[~, ind] = max(occurences);
largest_cluster = clusters(ind);

% Plot clustered coefficents 
[~, cluster_order] = sort(T);
ses_coef_sorted = ses_coef(cluster_order, :);
figure; imagesc(ses_coef_sorted'); colorbar;
xlabel('Sessions'); ylabel('Coefficients');
img_out = strcat(hc_distance(1), '_', hc_distance(2), ...
    '_coef_sorted.png');
saveas(gcf, fullfile(path_img_out, img_out));

% Plot objective function vs. # of clusters 
figure; plot(Z(:, 3));
ylabel('Objective Function'); xlabel('# of clusters');
img_out = strcat(hc_distance(1), '_', hc_distance(2), ...
    '_objective_func.png');
saveas(gcf, fullfile(path_img_out, img_out));

% Plot cluster assignment 
T_ses = reshape(T, [3 23]);
figure; imagesc(T_ses); colorbar
set(gcf, 'PaperUnits', 'inches');
x_width = 18; y_width = 2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
img_out = strcat(hc_distance(1), '_', hc_distance(2), ...
    '_cluster_assignment.png');
saveas(gcf, fullfile(path_img_out, img_out));

% Create a dendrogram plot of Z
% This is in case we wanted three clusters
% In reality, we will have to decide on 
% which # of clusters we want
cutoff = median([Z(end-2, n_clusters) Z(end-1, n_clusters)]);
figure; dendrogram(Z, n_sessions, 'ColorThreshold', cutoff);
img_out = strcat(hc_distance(1), '_', hc_distance(2), '_dendrogram.png');
saveas(gcf, fullfile(path_img_out, img_out));

% Select sessions belonging to the largest_cluster
sessions_in = (T == largest_cluster);
sessions_out = (T ~= largest_cluster);

% Compute and plot correlation between inliers  
corr_in = triu(corr(ses_coef(:, logical(sessions_in))), 2);
figure; imagesc(corr_in); colorbar;
msg = strcat('Avg. Corr:', " ", num2str(mean(corr_in(:))));
text(2, size(corr_in, 2), msg, 'FontSize', 14);
xlabel('Inlier Sessions', 'FontSize', 14);
ylabel('Inliear Sessions', 'FontSize', 14);
set(gcf, 'PaperUnits', 'inches');
x_width = 12; y_width = 11;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
img_out = strcat(hc_distance(1), '_', hc_distance(2), ...
    '_inlier_correlation.png');
saveas(gcf, fullfile(path_img_out, img_out));

% Compute and plot orrelation between outliers  
corr_out = triu(corr(ses_coef(:, logical(sessions_out))), 2);
figure; imagesc(corr_out); colorbar;
msg = strcat('Avg. Corr:', " ", num2str(mean(corr_out(:))));
text(2, size(corr_out, 2), msg, 'FontSize', 14);
xlabel('Outlier Sessions', 'FontSize', 14); 
ylabel('Outlier Sessions', 'FontSize', 14);
set(gcf, 'PaperUnits', 'inches');
x_width = 12; y_width = 11;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
img_out = strcat(hc_distance(1), '_', hc_distance(2), ...
    '_outlier_correlation.png');
saveas(gcf, fullfile(path_img_out, img_out));

end