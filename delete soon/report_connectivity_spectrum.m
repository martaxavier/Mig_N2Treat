%------------------------------------------------------------
% Connectivity spectrum 
%------------------------------------------------------------

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
    Dtf_f = signal; 

    my_title = strcat(upper(metric), ' between channels', ...
        " ", id_chans(plotting_channel), ' and', " ", ...
        id_chans(c));

    figure('Name', my_title);

    imagesc((1:size(signal,2)) ./ fs, ...
        1:length(f_vector), squeeze(Dtf_f));

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

