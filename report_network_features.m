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
        ' average (through time)', " ", upper(net_metric), ...
        ' for the', " ", id_bands(b), ' band');
    figure('Name', my_title);
    topoplot(signal, chanlocs, topo_settings{:});
    title(my_title);
    colorbar; caxis([-max_signal_abs max_signal_abs]); 
    
    img_out = strcat(upper(full_metric), ...
        '_', id_bands(b), '_TOPO_Avg.png');
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
        num2str(delays(plotting_delay)), 'sec', ...
        id_chans(plotting_channel),'.png');
    saveas(gcf,fullfile(path_img_out(s), img_out));

end % looping through bands

        