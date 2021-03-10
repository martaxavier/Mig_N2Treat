%---------------------------------------------------------    
% Plot random connectome 
%---------------------------------------------------------

% Representation of the connectome at a random time-point, 
% for each frequency band 

% Choose a time-point randomly 
rand_pnt = floor(1 + rand(1)*(n_pnts - 1));

% Types of connectivity/p-value
% matrices to plot 
my_files = {'Connectome', ...
            'PValues', ...
            'ConnectomeSig', ...
            'ConnectomeTopo'};

my_titles = {'Connectome', ...
            'P-values', ...
            'Connectome (Stat. Filt.)', ...
            'Connectome (Topo. Filt.)'};
        
% Go through bands         
for b = 1 : n_bands
    
    % Plot the connectome before and after statistical filtering, 
    % and after topological filtering; plot the matrix of p-values
    my_signals = zeros([size(squeeze(conspec(rand_pnt, :, :, b))) 4]);
    my_signals(:, :, 1) = squeeze(conspec(rand_pnt, :, :, b));
    my_signals(:, :, 2) = squeeze(p_values(rand_pnt, :, :, b));
    my_signals(:, :, 3) = squeeze(conspec_sig(rand_pnt, :, :, b));
    my_signals(:, :, 4) = squeeze(conspec_topo(rand_pnt, :, :, b));    
      
    % Go through types of signals 
    for t = 1 : size(my_signals, 3)

        % Generate current connectivity/p-value matrix 
        my_title = strcat(my_titles(t), ' at time-point', " ", ...
            num2str(rand_pnt),' , for the', " ", id_bands(b), ' band');
        fig = figure('Name', char(my_title));
        fig.Position(3 : 4) = fig.Position(3 : 4)*5;
        imagesc(squeeze(my_signals(:, :, t))); colorbar; 
        xlabel('Channels','FontSize', 24); 
        ylabel('Channels','FontSize', 25);
        xticks(1 : n_chans); yticks(1 : n_chans);
        xticklabels(cellstr(id_chans));
        yticklabels(cellstr(id_chans));
        set(gcf, 'PaperUnits', 'inches');
        x_width = 12; y_width = 11;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

        % Save generated images in output path 
        img_out = strcat(upper(metric), ...
            '_', my_files(t), '_', id_bands(b),'.png');
        saveas(gcf, char(fullfile(path_img_out(s), img_out)));
        
    end

end

%---------------------------------------------------------    
% Plot average connectome 
%---------------------------------------------------------
  
% Representation of the average connectome across time and 
% frequency bands 

% Plot the connectome before and after statistical filtering, 
% and after topological filtering; plot the matrix of p-values
my_signals = zeros([size(squeeze(mean(squeeze(mean(conspec, 1)), 3))) 4]);
my_signals(:, :, 1) = squeeze(mean(squeeze(mean(conspec, 1)), 3));
my_signals(:, :, 2) = squeeze(mean(squeeze(mean(p_values, 1)), 3));
my_signals(:, :, 3) = squeeze(mean(squeeze(mean(conspec_sig, 1)), 3));
my_signals(:, :, 4) = squeeze(mean(squeeze(mean(conspec_topo, 1)), 3));

my_files = {'Connectome_Avg', ...
            'PValues_Avg', ...
            'ConnectomeSig_Avg', ...
            'ConnectomeTopo_Avg'};
   
% Go through types of signals 
for t = 1 : size(my_signals, 3)

    % Generate current connectivity/p-value matrix 
    my_title = strcat('Average', " ", my_titles(t), ...
        ' across time-points and frequency bands');
    fig = figure('Name', my_title);
    fig.Position(3 : 4) = fig.Position(3 : 4)*5;
    imagesc(squeeze(my_signals(:, :, t))); colorbar; 
    xlabel('Channels','FontSize', 24); 
    ylabel('Channels','FontSize', 24);
    xticks(1 : n_chans); yticks(1 : n_chans);
    xticklabels(cellstr(id_chans));
    yticklabels(cellstr(id_chans));
    set(gcf, 'PaperUnits', 'inches');
    x_width = 12; y_width = 11;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

    % Save generated images in output path 
    img_out = strcat(upper(metric), ...
        '_', my_files(t), '.png');
    saveas(gcf, char(fullfile(path_img_out(s), img_out)));

end
