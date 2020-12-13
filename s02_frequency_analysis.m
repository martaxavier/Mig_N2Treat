import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 
 
% ------------------------------------------------------------
% Power Spectral Density Topographies 
% ------------------------------------------------------------ 

% Create a title for the report
my_title = 'POWER TOPOGRAPHIES';
H1 = get_report_heading(1,my_title);
add(R,H1)

%---------------------------------------------------------    
% Go through metrics
%---------------------------------------------------------    

for m = 1 : length(metrics)

    metric = metrics(m);
    get_metric_pars;
        
    metric_data_in = strcat...
        (eeg_metric,'_', data_in);
    
    my_title = upper(metric);
    H2 = get_report_heading(2,my_title);
    add(R,H2);

    %-------------------------------------------------------
    % Go through subjects 
    %-------------------------------------------------------

    for s = 1 : length(subjects)
        
        % Define current subject 
        subject = subjects(s);
        
        disp(strcat('Performing frequency analysis for', ... 
            " ", subject,', metric'," ", metric, ' ...'));
        
        my_title = subject;
        H3 = get_report_heading(3,my_title);
        add(R,H3);

        % Create subject ouput directory 
        % if not existent 
        if ~exist(path_img_out(s), 'dir')
            mkdir(path_img_out(s))
        end

        %-----------------------------------------------------
        % Load data
        %-----------------------------------------------------

        features = dlmread(fullfile(path_data_in(s), ...
            metric_data_in));
        n_pnts = size(features,1);
        features = reshape(features,[n_pnts,n_chans,n_bands]);

        %-----------------------------------------------------
        % Plot time-averaged topographies 
        %-----------------------------------------------------

        % These are power spectral density values (Fourier
        % transform squared) so they are all positive, but 
        % less than 1. Hence, they log values will be negative. 
        % Higher (or less negative) values correspond to higher
        % power spectral density, and vice-versa.
        % Log values are used to increase sensitivity in power 
        % differences between lower power values. This is because 
        % the color gradient of the topographic images will vary 
        % linearly with the log of the power. 
        features_log = squeeze(mean(10.*log10(features)));
        if length(size(features)) == 2
            features_log = features_log';
        end

        subp_cols = ceil(n_bands/2);
        subp_rows = ceil(n_bands/subp_cols);
        fig = figure();

        for b = 1 : n_bands

            % Minimum and maximum intensity values 
            min_val = min(features_log(:,b));
            min_int = min_val - 0.03*abs(min_val);
            max_val = max(features_log(:,b));
            max_int = max_val + 0.03*abs(max_val);

            subplot(subp_rows,subp_cols,b);
            title(id_bands(b)); 
            topoplot(features_log(:,b),chanlocs,'electrodes',...
               'labels','whitebk','on','gridscale',500);
            colorbar; caxis([min_int max_int]);       

        end          

        %Increase your figure pixel resolution
        fig.Position(3:4) = fig.Position(3:4)*5;

        img_out = strcat(upper(metric),'_LOGTOPOPOWER.png');
        source = fullfile(path_img_out(s),img_out);
        saveas(fig,source); I = Image(source);
        I.Style={ScaleToFit(true),HAlign('center')};
        add(R,I);

    end % finish looping through subjects

end % finish looping through metrics 

% ------------------------------------------------------------
% Spectral Power Time-courses  
% ------------------------------------------------------------ 

% Create a title for the report
my_title = 'POWER TIME-COURSES';
H1 = get_report_heading(1,my_title);
add(R,H1);

%---------------------------------------------------------    
% Go through metrics
%---------------------------------------------------------    

for m = 1 : length(metrics)

    metric = metrics(m);
    get_metric_pars;
    
    my_title = upper(metric);
    H2 = get_report_heading(2,my_title);
    add(R,H2);

    %-------------------------------------------------------
    % Go through subjects 
    %-------------------------------------------------------

    for s = 1 : length(subjects)

        % Define current subject 
        subject = subjects(s);

        my_title = subject;
        H3 = get_report_heading(3,my_title);
        add(R,H3);
        
        report_channels = [4 10];

        %-----------------------------------------------------
        % Report TF time-courses 
        %-----------------------------------------------------

        for c = 1 : length(report_channels)

            chan = report_channels(c);

            my_tytle = upper(strcat(id_bands(1), ...
                ', channel'," ",num2str(chan)));
            H4 = get_report_heading(4,my_title);
            add(R,H4);

            img_out = strcat(id_bands(1),'5sec', ...
                num2str(chan),'chan.png');
            source = fullfile(path_img_in(s),img_out);
            I = Image(source);
            I.Style={ScaleToFit(true),HAlign('center')};
            add(R,I);   

        end  % finish looping through channels        

    end % finish looping through subjects

end % finish looping through metrics 
