%------------------------------------------------------------    
% Go through subjects
%------------------------------------------------------------
for s = 1 : length(subjects)

    % Define current subject 
    subject = subjects(s);

    disp(strcat('Extracting markers for'," ", subject,' ...'));

    %---------------------------------------------------------    
    % Load data 
    %--------------------------------------------------------- 

    % Create output directories if non-existent 
    if ~exist(path_data_out(s), 'dir'); mkdir(path_data_out(s)); end  

    % Load eeg dataset of specified subject 
    load(fullfile(path_data_in(s),data_in));
    
    % Extract data from dataset
    data = EEG.data;                     % eeg data matrix 
    n_pnts = size(data,2);               % number of datapoints

    % Extract event information 
    event = struct2cell(EEG.event);
    type = event(1,:,:); type = squeeze(type);
    latency = event(2,:,:); latency=squeeze(latency);
    latency = cell2mat(latency);

    % Extract channel information 
    chanlocs = struct2cell(EEG.chanlocs);
    labels = chanlocs(1,:,:); labels = squeeze(labels);

    % Find EEG sample that corresponds to the beginning 
    % (first_eeg) and end (last_eeg) of the simultaneous 
    % BOLD acquisition 
    idxs_scan = find(ismember(type,markers_task));
    first_eeg = round(latency(idxs_scan(1))); 
    last_eeg = round(latency(idxs_scan(end)));
    
     if sub_task ~= ""
         
        % Find EEG samples that correspond to the beginning 
        % and end of the epochs of the current sub_task 
        idxs_sub_task_start = find(ismember(type,markers_sub_task_start));
        idxs_sub_task_stop =  find(ismember(type,markers_sub_task_stop));
        sub_task_start_eeg = round(latency(idxs_sub_task_start));
        sub_task_stop_eeg = round(latency(idxs_sub_task_stop));

        n_sub_task_start_eeg = length(sub_task_start_eeg);
        n_sub_task_stop_eeg = length(sub_task_stop_eeg);

        % Defiene a design matrix that is one for the samples that belong 
        % to the current sub_task, and zero for the samples that don't
        design_table = cat(2,sub_task_start_eeg,ones(n_sub_task_start_eeg,1));
        design_table = cat(1,design_table,cat(2,sub_task_stop_eeg, ...
            zeros(n_sub_task_stop_eeg,1)));

        [~,idx] = sort(design_table(:,1));
        design_table = design_table(idx,:);
        n_markers = size(design_table,1);

        sub_task_design = zeros(1,n_pnts);

        for m = 1 : n_markers

            marker_start = design_table(m,1);

            if m == n_markers; marker_stop = n_pnts;
            else;marker_stop = design_table(m+1,1)-1;end

            marker_value = design_table(m,2);
            sub_task_design(marker_start:marker_stop) = marker_value;

        end
        
     end

end
