import mlreportgen.dom.*;
import mlreportgen.report.*;
iptgetpref('ImshowInitialMagnification'); 

% Create a title for the report
if flag.report ~= 0
    my_title = 'POWER RELIABILITY ANALYSIS';
    H1 = get_report_heading(1, my_title);
    add(R, H1);
end
 
n_metrics = length(metrics);
n_sessions = length(sessions);
n_subjects = length(subjects);

%---------------------------------------------------------    
% Retreive subjects average power 
%---------------------------------------------------------    

subj_power = cell(n_metrics, 1);                         

% Go through metrics 
for m = 1 : length(metrics)

    metric = metrics(m);
    get_metric_pars;
    
    if ~contains(metric, power_metrics)
        continue;
    end    
        
    % Pre-allocate power matrix 
    power_all = zeros(n_chans*n_bands, n_subjects, n_sessions);

    % Go through sessions 
    for se = 1 : n_sessions
            
        % Go through subjects 
        for s = 1 : length(subjects)

            % Load subject power data
            load(fullfile(path_data_in(s, se), ...
                strcat(metric, '_', data_in)));

            % Save average power for current session and subject
            power_all(:, s, se) = power(:);
            
        end % subjects 
        
    end % sessions
    
    subj_power{m, 1} = reshape(power_all, [n_chans*n_bands n_subjects*n_sessions]);

end % metrics 

% Go through metrics 
for m = 1 : n_metrics

    metric = metrics(m);
    get_metric_pars;
    n_features = prod(dim);  
        
    % Broadcast the current pipeline stage 
    disp(strcat('Performing power reliability', ...
        ' analysis for metric'," ", metric, ' ...'));
        
    % Add metric heading 
    % to the report 
    if flag.report ~=0
        
        my_title = upper(metric);
        H2 = get_report_heading(2, my_title);
        add(R, H2);
        
    end   
    
    % ------------------------------------------------------
    % Reliability tests tests 
    % ------------------------------------------------------    
    
    % Compute measures of intra- and inter-subject variability,
    % as well as the interclass correlation coefficient for each 
    % EEG feature 
    reliab = reliability_test(cell2mat(subj_power(m, ...
        1)), n_subjects, n_sessions, reliab_metric, ...
        fullfile(path_data_out, strcat(metric, '_', data_out)));
    
    half_corr = reliability_split_half_test(cell2mat(subj_power(m, 1)), ...
        strcat(metric, '_deconv'), n_subjects, n_sessions, ...
        n_iter_reliab, path_pars, path_data_out, path_img_out);    
    
    % ------------------------------------------------------
    % Plot reliability results   
    % ------------------------------------------------------    
    
    if flag.report ~= 0
        
        % Add report title 
        my_title = upper(reliab_metric);
        H3 = get_report_heading(3, my_title);
        add(R, H3);

        report_reliability(cell2mat(subj_power(m, 1)), ...
            reliab, strcat(metric, '_deconv'), R, path_data_out, ...
            path_img_out, path_pars); 
        
    end 
    
end % metrics

