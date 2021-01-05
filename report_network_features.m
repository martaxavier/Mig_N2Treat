%---------------------------------------------------------    
% Plot random network
%---------------------------------------------------------

% NOTE: conspec_sig is n_pnts x n_chans x n_chans x n_bands 

% Representation of the network at a random time point
% For each frequency band

rand_pnt = floor(1 + rand(1)*(n_pnts - 1));

for b = 1 : n_bands
    
    conspec_topo_rand = squeeze(conspec_topo(rand_pnt, :, :, b)); 
    
    my_title = strcat('Connectome at a random', ...
        ' time-point, for the', id_bands(b), ' band');
    fig = figure('Name', my_title);
    fig.Position(3:4) = fig.Position(3:4)*5;
    imagesc(conspec_topo_rand); colorbar; 
    xlabel('Channels','FontSize',24); 
    ylabel('Channels','FontSize',24);
    xticks(1:n_chans); yticks(1:n_chans);
    xticklabels(cellstr(id_chans));
    yticklabels(cellstr(id_chans));
    set(gcf, 'PaperUnits', 'inches');
    x_width = 12; y_width = 11;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
 
    img_out = strcat(upper(con_metric), ...
        '_Connectome_',id_bands(b),'.png');
    saveas(gcf,fullfile(path_img_out(s), img_out));

end

%---------------------------------------------------------    
% Plot average networks 
%---------------------------------------------------------

% Representation of the average network across time and 
% frequency bands 
conspec_sig_avg = squeeze(mean(squeeze(mean(conspec_sig, 1)), 3));

my_title = strcat('Average connectome,', ...
    ' throughout bands and time, filtered');
fig = figure('Name', my_title);
fig.Position(3:4) = fig.Position(3:4)*5;
imagesc(conspec_sig_avg); colorbar;
xlabel('Channels','FontSize',20); 
ylabel('Channels','FontSize',20);
xticks(1:n_chans); yticks(1:n_chans);
xticklabels(cellstr(id_chans));
yticklabels(cellstr(id_chans));
set(gcf, 'PaperUnits', 'inches');
x_width = 12; y_width = 11;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
    
img_out = strcat(upper(con_metric), ...
    '_Connectome_avg.png');
saveas(gcf,fullfile(path_img_out(s), img_out));

%---------------------------------------------------------    
% Plot average network measures topographies 
%---------------------------------------------------------

% Average topographies of the network measures
% across time 

topo_settings = {'electrodes', 'labels', ...
    'whitebk', 'on', 'gridscale', 300};
    
for b = 1 : n_bands
    
    signal = squeeze(eeg_features_norm(:, :, b));
    signal = squeeze(mean(signal, 1));
    
    max_signal_abs = max(max(signal), ...
    abs(min(signal)));
    
    my_title = strcat('Topographic map of', ...
        'average (through time)', " ", upper(net_metric), ...
        ' for the', " ", id_bands(b), ' band');
    topoplot(signal, chanlocs, topo_settings{:});
    colorbar; caxis([-max_signal_abs max_signal_abs]); 
    
    img_out = strcat(upper(full_metric), ...
        '_', id_bands(b), '_topography_avg.png');
    saveas(gcf,fullfile(path_img_out(s), img_out));
    
end

%---------------------------------------------------------    
% Plot all feature (network measure) time-series  
%---------------------------------------------------------
        
for b = 1 : n_bands

    signal = eeg_features_norm(:, plotting_channel, b);
    signal_delayed = eeg_features_delayed_norm...
        (:, plotting_channel, plotting_delay, b);

    my_title = strcat(id_bands(b), " ", upper(con_metric), ...
        " ", upper(net_metric), ' of channel', " ", ...
        id_chans(plotting_channel), ', delay', " ", ...
        num2str(delays(plotting_delay)), 's');

    figure('Name',my_title)

    plot(time, signal); hold on; 
    plot(time, signal_delayed);

    title(my_title); 
    xlabel('Time(s)'); 
    ylabel(upper(strcat(con_metric, ...
        '{_', net_metric, '}')));
    legend(strcat(id_bands(b),' feature'), ...
        strcat(id_bands(b),' feature,', ...
        plotting_shift,', delay'," " , ...
        num2str(delays(plotting_delay)), 's'));

    img_out = strcat(upper(full_metric), '_', id_bands(b), ...
        num2str(delays(plotting_delay)), 'sec_', ...
        id_chans(plotting_channel),'.png');
    saveas(gcf,fullfile(path_img_out(s), img_out));

end % looping through bands

        