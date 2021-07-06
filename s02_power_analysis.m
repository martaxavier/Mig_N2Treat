import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
 
% ------------------------------------------------------------
% Power Spectral Density Topographies 
% ------------------------------------------------------------ 


% Create a title for the report
if flag.report ~=0 
    % Create a title for the report
    my_title = strcat('POWER ANALYSIS -', ...
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
        add(R,H2);    
    end

    %-------------------------------------------------------
    % Go through subjects 
    %-------------------------------------------------------

    for s = 1 : length(subjects)
        
        % Define current subject 
        subject = subjects(s);
        
        get_metric_pars
        
        disp(strcat('Performing power analysis for', ... 
            " ", subject, ', metric', " ", metric, ' ...'));
        
        if flag.report ~=0         
            my_title = subject;
            H3 = get_report_heading(3, my_title);
            add(R, H3);         
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
        % Estimate Power Spectral Density w/ Welch's method 
        %----------------------------------------------------- 

        f_vector = linspace(f_min, f_max, n_freq);
        [pxx, ~] = pwelch(data', [], [], f_vector, fs_eeg);
        power = average_frequency(pxx, f_vector, bands, 1);
        
        % These are power spectral density values (Fourier
        % transform squared) so they are all positive, but 
        % less than 1. Hence, they log values will be negative. 
        % Higher (or less negative) values correspond to higher
        % power spectral density, and vice-versa.
        % Log values are used to increase sensitivity in power 
        % differences between lower power values. This is because 
        % the color gradient of the topographic images will vary 
        % linearly with the log of the power.        
        power_log = 10.*log10(power);

        if size(power, 1) > 1
            power = power'; power_log = power_log';
        end
        
        % Save power 
        save(fullfile(path_data_out(s, se), ...
            strcat(metric, '_', data_out)), 'power');
        
        % Continue to skip plots
        if flag.report == 0; continue; end
        
        %-----------------------------------------------------
        % Plot Power topographies 
        %-----------------------------------------------------

        subp_cols = ceil(n_bands/2);
        subp_rows = ceil(n_bands/subp_cols);
        fig = figure();
            
        for b = 1 : n_bands
             
            % Minimum and maximum intensity values 
            min_val = min(power_log(:, b)) - 0.1*min(power_log(:, b));
            max_val = max(power_log(:, b)) + 0.1*max(power_log(:, b));       

            subplot(subp_rows, subp_cols, b);
            title(id_bands(b)); 
            topoplot(power_log(:, b), chanlocs, 'electrodes', ...
               'labels', 'whitebk', 'on', 'gridscale', 100);
            colorbar; caxis([min_val max_val]);       

        end          

        %Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;

        img_out = strcat(upper(metric), '_LOGTOPOPOWER.png');
        source = fullfile(path_img_out(s, se), img_out);
        saveas(fig, source); I = Image(source);
        I.Style = {ScaleToFit(true), HAlign('center')};
        add(R,I);

    end % finish looping through subjects

end % finish looping through metrics 

% Return to skip plots 
if flag.report == 0; return; end

% ------------------------------------------------------------
% Plot Spectral Power Time-courses  
% ------------------------------------------------------------ 

% Create a title for the report
my_title = 'POWER TIME-COURSES';
H1 = get_report_heading(1, my_title);
add(R, H1);

%---------------------------------------------------------    
% Go through metrics
%---------------------------------------------------------    

for m = 1 : length(metrics)

    metric = metrics(m);
    get_metric_pars;
    
    my_title = upper(metric);
    H2 = get_report_heading(2, my_title);
    add(R, H2);

    %-------------------------------------------------------
    % Go through subjects 
    %-------------------------------------------------------

    for s = 1 : length(subjects)

        % Define current subject 
        subject = subjects(s);

        get_metric_pars; 
        
        my_title = subject;
        H3 = get_report_heading(3, my_title);
        add(R,H3);
        
        report_channels = [4 10];

        %-----------------------------------------------------
        % Report TF time-courses 
        %-----------------------------------------------------

        for c = 1 : length(report_channels)

            chan = report_channels(c);

            my_tytle = upper(strcat(id_bands(1), ...
                ', channel', " ", num2str(chan)));
            H4 = get_report_heading(4, my_title);
            add(R, H4);

            img_out = strcat(id_bands(1), '5sec', ...
                num2str(chan), 'chan.png');
            source = fullfile(path_img_in(s, se), img_out);
            I = Image(source);
            I.Style={ScaleToFit(true), HAlign('center')};
            add(R,I);   

        end  % finish looping through channels        

    end % finish looping through subjects

end % finish looping through metrics 
