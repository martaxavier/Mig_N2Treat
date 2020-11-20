%------------------------------------------------------------
% Connectivity spectrum 
%------------------------------------------------------------

alpha = (f_max/f_min)^(1/n_freq) - 1;

% Full connectivity 
if size(squeeze(dim)) == 4
  
    % Create symmetrical matrix 
    conspec_T = permute(conspec,[1,2,4,3]);
    conspec_sym = conspec_T + conspec;
    
    % Plot connectivity spectrum (not normalized), between the specified 
    % plotting channel and all other channels, throughout time     
    for c = 1 : n_chans
        
        if c == plotting_channel
            continue
        end

        % Define current signal from the connectivity spectrum 
        signal = squeeze(conspec_sym(:, :, plotting_channel, c));

        % Define colorscale 
        %Dtf_f = log(abs(signal) + .005); 
        Dtf_f = signal; 
        
        my_title = strcat(upper(metric), ' between channels', ...
            " ", id_chans(plotting_channel), ' and', " ", ...
            id_chans(c));

        figure('Name', my_title);

        imagesc((1:size(signal,2)) ./ fs, ...
            1:length(f_vector), squeeze(Dtf_f));

        Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) ...
            - log(f_vector(1))) ./ log(1 + alpha) + 2;
        hold on; plot([1 n_pnts],[Fplot', Fplot'],'k');
        hold off

        title(my_title);
        set(gca,'YTick',Fplot);
        set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90],'FontSize',12)
        ylabel('Frequency (Hz)','FontSize', 26);
        xlabel('Time (s)','FontSize',26);
        colorbar;
        
        img_out = strcat(upper(metric), '_', id_chans...
            (plotting_channel), 'vs', id_chans(c), '.png');
        saveas(gcf,fullfile(path_img_out(s), img_out));

    end % finish looping through channels 
 
% Graph metrics     
else
   
    % Plot connectivity spectrum (not normalized),of the 
    % specified plotting channel, throughout time     
    for c = 1 : n_chans

        % Define current signal from the connectivity spectrum 
        signal = squeeze(conspec(:, :, plotting_channel));

        % Define colorscale 
        %Dtf_f = log(abs(signal) + .001); 
        Dtf_f = signal;
        
        C_metric = extractBefore(metric,'_');
        G_metric = extractAfter(metric,'_');

        my_title = strcat(upper(C_metric), " ", upper(G_metric), ...
            'between channels', " ", id_chans(plotting_channel), ...
            ' and', " ", id_chans(c));

        figure('Name', my_title);

        imagesc((1:size(signal,2)) ./ fs, ...
            1:length(f_vector), squeeze(Dtf_f));

        Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) ...
            - log(f_vector(1))) ./ log(1 + alpha) + 2;
        hold on; plot([1 n_pnts],[Fplot', Fplot'],'k');
        hold off

        title(my_title);
        set(gca,'YTick',Fplot);
        set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90],'FontSize',12)
        ylabel('Frequency (Hz)','FontSize', 26);
        xlabel('Time (s)','FontSize',26);
        colorbar;
        
        img_out = strcat(upper(metric), '_', id_chans...
            (plotting_channel), 'vs', id_chans(c), '.png');
        saveas(gcf,fullfile(path_img_out(s), img_out));
        
    end % finish looping through channels 
    
end
    
%------------------------------------------------------------
% Features (normalized, convolved)
%------------------------------------------------------------

% Plot features (exactly as in compute_power_features) 
% but run through all channels in case of icoh 
if ~isempty(eeg_shift) 

    % Full connectivity 
    if size(squeeze(dim)) == 4
    
        % Create symmetrical matrices  
        features_T = permute(eeg_features_norm,[1,3,2,4]);
        features_delayed_T = permute(eeg_features_delayed_norm, ...
            [1,3,2,4,5]); 
        
        features_sym = features_T + eeg_features_norm;
        features_delayed_sym = features_delayed_T + ...
            eeg_features_delayed_norm;
        
        for c = 1 : n_chans
            
            if c == plotting_channel
                return
            end

            for b = 1 : n_bands

                signal = features_sym(:, plotting_channel, c, b);
                signal_delayed = features_delayed_sym...
                    (:, plotting_channel, c, plotting_delay, b);

                my_title = strcat(id_bands(b), " ", upper(metric), ...
                    ' between channel', " ", id_chans(plotting_channel), ...
                    ' and channel', " ", id_chans(c), ' delay', " ", ...
                    num2str(delays(plotting_delay)), 's');

                figure('Name',my_title)

                plot(time, signal); hold on; 
                plot(time, signal_delayed);

                title(my_title); 
                xlabel('Time(s)'); ylabel(upper(metric));
                legend(strcat(id_bands(b),' feature'), ...
                    strcat(id_bands(b),' feature,', ...
                    plotting_shift,', delay'," " , ...
                    num2str(delays(plotting_delay)), 's'));

                img_out = strcat(upper(metric),'_',id_bands(b), ...
                    num2str(delays(plotting_delay)),'sec_', ...
                    id_chans(plotting_channel),'vs',id_chans(c), ...
                    '.png');
                saveas(gcf,fullfile(path_img_out(s),img_out));

            end % looping through bands

        end % looping through channels
        
    % Graph metrics     
    else
       
        for c = 1 : n_chans     

            for b = 1 : n_bands

                signal = features_sym(:, plotting_channel, b);
                signal_delayed = features_delayed_sym...
                    (:, plotting_channel, plotting_delay, b);

                my_title = strcat(id_bands(b), " ", upper(C_metric), ...
                    " ", upper(G_metric), ' between channel', " ", ...
                    id_chans(plotting_channel), ' and channel', " ", ...
                    id_chans(c), ' delay', " ", num2str(delays...
                    (plotting_delay)), 's');

                figure('Name',my_title)

                plot(time, signal); hold on; 
                plot(time, signal_delayed);

                title(my_title); 
                xlabel('Time(s)'); ylabel(upper(metric));
                legend(strcat(id_bands(b),' feature'), ...
                    strcat(id_bands(b),' feature,', ...
                    plotting_shift,', delay'," " , ...
                    num2str(delays(plotting_delay)), 's'));

                img_out = strcat(upper(metric), '_', id_bands(b), ...
                    num2str(delays(plotting_delay)), 'sec_', ...
                    id_chans(plotting_channel), 'vs', id_chans(c), ...
                    '.png');
                saveas(gcf,fullfile(path_img_out(s), img_out));

            end % looping through bands

        end % looping through channels
        
    end

end
        
if size(squeeze(dim)) == 4
    return
end

% Plot topoplots of weighted node degree, at each frequency 
% band averaged through time 

%  
 