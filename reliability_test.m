function [reliab] = reliability_test(data, n_subjects, n_sessions, reliab_metric, path_data_out)

%-----------------------------------------------------------
% Perform variability tests 
%-----------------------------------------------------------

n_features = size(data, 1);
data = reshape(data, [n_features, ...
    n_subjects, n_sessions]);

% Allocate 
reliab = zeros(n_features, 1);

for f = 1 : n_features

    % Data of current feature
    data_f = squeeze(data(f, :, :)); 
    
    % If the formula contains the absolute value of the mean,
    % the measure is called relative standard deviation instead
    % CV higher than 100 just means that the standard 
    % deviation is higher than the mean 
    % Cv negative just means that the mean is negative 
    
    switch reliab_metric
        
        case 'cv_inter'
            
            % CV-inter is computed as the average cv-inter across sessions
            reliab(f) = mean(100*(std(data_f) ./  mean(data_f)));
            
        case 'cv_intra_1'
                     
            % CV-intra (version 1)
            mean_subj = mean(data_f, 2);    
            std_subj = sqrt(mean(std(data_f, 1, 2)));
            reliab(f) = 100*(std_subj / mean(mean_subj));
            
        case 'cv_intra_2'
            
            % CV-intra (version 2)
            mean_subj = mean(data_f, 2);  
            cv_subj = mean(100*(std(data_f, 1, 2) ./ mean_subj));
            reliab(f) = sqrt(mean(cv_subj.^2));
    
        case 'icc'        
               
            % ICC 
            ss_total = var(data_f(:)) *(n_subjects*n_sessions - 1);
            MSR = var(mean(data_f, 2)) * n_sessions;
            MSC = var(mean(data_f, 1)) * n_subjects;
            MSE = (ss_total - MSR *(n_subjects - 1) - MSC * ...
                (n_sessions - 1))/ ((n_subjects - 1) * (n_sessions - 1));
            reliab(f) = (MSR - MSE) / (MSR + (n_sessions-1)*MSE);   

    end
    
end % features 

% Negative ICC results will be set to zero for 
% since in this case negative values are not 
% relevant for the analysis
if strcmp(reliab_metric, 'icc')
    reliab(reliab<0) = 0; 
end

% Save results 
save(path_data_out, 'reliab');

end