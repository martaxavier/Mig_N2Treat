import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
 
% ------------------------------------------------------------
% Power Spectral Density Topographies 
% ------------------------------------------------------------ 

% Create a title for the report
if flag.report ~=0 
    my_title = strcat('CONNECTIVITY ANALYSIS -', ...
        " ", upper(session));
    H1 = get_report_heading(1, my_title);
    add(R, H1);
end

%---------------------------------------------------------    
% Go through metrics
%---------------------------------------------------------    

for m = 1 : length(metrics)

    metric = metrics(m);
    get_metric_pars;
    
    if flag.report ~=0 
        my_title = upper(metric);
        H2 = get_report_heading(2, my_title);
        add(R, H2);
    end

    %-------------------------------------------------------
    % Go through subjects 
    %-------------------------------------------------------

    for s = 1 : length(subjects)
        
        % Define current subject 
        subject = subjects(s);
        
        get_metric_pars
        
        disp(strcat('Performing connectivity analysis for', ... 
            " ", subject, ', metric', " ", metric, ' ...'));
        
        if flag.report ~=0 
            my_title = subject;
            H3 = get_report_heading(3, my_title);
            add(R,H3);
        end

        % Create subject ouput directories 
        % if not existent 
        if ~exist(path_data_out(s, se), 'dir')
            mkdir(path_data_out(s, se))
        end
        
        if ~exist(path_img_out(s, se), 'dir')
            mkdir(path_img_out(s, se))
        end

        %-----------------------------------------------------
        % Load data
        %-----------------------------------------------------

        % Load eeg dataset of specified subject 
        load(fullfile(path_data_in(s, se), data_in));

        first_last = dlmread(fullfile(path_markers_in(s, se), markers_in));
        
        % Assign first and last 
        % EEG samples
        first_eeg = first_last(1); 
        last_eeg = first_last(end); 

        % Extract EEG data from dataset
        data = EEG.data;   
        data = data(:, first_eeg : last_eeg);
        
        %-----------------------------------------------------
        % Estimate static functional connectivity 
        %----------------------------------------------------- 

        f_vector = linspace(f_min, f_max, n_freq);
        max_freq_res = max(diff(f_vector)); 
        
        [cxy, p_values, ~] = get_connectivity(data, fs_eeg, ...
            max_freq_res, con_metric);
        [decision, p_thresh_corrected, ~, ~] = fdr_bh(p_values, ...
            thresh_fdr, 'pdep', 'no');
        cxy(~decision) = 0;
        cxy(isnan(cxy)) = 0;

        con = average_frequency(cxy, f_vector, bands, 2);
        
        % Save connectivity 
        save(fullfile(path_data_out(s, se), ...
            strcat(con_metric, '_', data_out)), 'con');
        
        % Convert into symmetric matrix 
        con = vec2tri(con, 1, n_chans, 'upper');
        con = tri2symmetric(con, 1, 2);        
        
        %-----------------------------------------------------
        % Estimate network metric  
        %-----------------------------------------------------
        
        con_sing = zeros([1 size(con)]);
        con_sing(1, :, :, :) = con; con_sing(isnan(con_sing)) = 0;
        Gxy = compute_network_metric(con_sing, net_metric);
        
        % Save network measure  
        net = Gxy;
        save(fullfile(path_data_out(s, se), ...
            strcat(metric, '_', data_out)), 'net');
        
        % Continue to skip plots 
        if flag.report == 0; continue; end
        
        %-----------------------------------------------------
        % Plot static connectivity matrices 
        %-----------------------------------------------------
        
        subp_cols = ceil(n_bands/2);
        subp_rows = ceil(n_bands/subp_cols);
        fig = figure();
        set(gcf, 'PaperUnits', 'inches');
        x_width = 12; y_width = 11;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
        
        for b = 1 : n_bands

            subplot(subp_rows, subp_cols, b);
            title(id_bands(b)); 
            imagesc(con(:, :, b)); colorbar; 
            xlabel('Channels','FontSize', 24); 
            ylabel('Channels','FontSize', 24);      

        end          

        % Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;

        img_out = strcat(upper(con_metric), '_CONNECTIVITY.png');
        source = fullfile(path_img_out(s, se), img_out);
        saveas(fig, source); I = Image(source);
        I.Style={ScaleToFit(true), HAlign('center')};
        add(R, I);        
        
        %-----------------------------------------------------
        % Plot network metric topographies  
        %-----------------------------------------------------
        
        subp_cols = ceil(n_bands/2);
        subp_rows = ceil(n_bands/subp_cols);
        fig = figure();
            
        for b = 1 : n_bands
             
            % Minimum and maximum intensity values 
            min_val = min(Gxy(:, b)) - 0.1*min(Gxy(:, b));
            max_val = max(Gxy(:, b)) + 0.1*max(Gxy(:, b));       

            subplot(subp_rows, subp_cols, b);
            title(id_bands(b)); 
            topoplot(Gxy(:, b), chanlocs, 'electrodes', ...
               'labels', 'whitebk', 'on', 'gridscale', 100);
            colorbar; caxis([min_val max_val]);       

        end  
        
        % Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;

        img_out = strcat(upper(net_metric), '_NETWORK.png');
        source = fullfile(path_img_out(s, se), img_out);
        saveas(fig, source); I = Image(source);
        I.Style={ScaleToFit(true), HAlign('center')};
        add(R, I);        

    end % finish looping through subjects

end % finish looping through metrics 