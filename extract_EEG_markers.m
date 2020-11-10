function [first_eeg,last_eeg,sub_task_design] = ...
    extract_EEG_markers(EEG,sub_task,markers_task, ...
    markers_sub_task_start,markers_sub_task_stop, sub_task_order)


% Extract data from dataset
data = EEG.data;                     % eeg data matrix 
n_pnts = size(data,2);               % number of datapoints

% Extract event information 
event = struct2cell(EEG.event);
type = event(1,:,:); type = squeeze(type);
latency = event(2,:,:); latency=squeeze(latency);
latency = cell2mat(latency);

% Find EEG sample that corresponds to the beginning 
% (first_eeg) and end (last_eeg) of the simultaneous 
% BOLD acquisition 
idxs_scan = find(ismember(type,markers_task));
first_eeg = round(latency(idxs_scan(1))); 
last_eeg = round(latency(idxs_scan(end)));

 if sub_task ~= ""
     
    % Find EEG samples that correspond to the  
    % beginning and end of the current task 
    task_start_eeg = round(latency...
        (ismember(type,markers_task_start)));
    task_stop_eeg = round(latency...
        (ismember(type,markers_task_stop)));
     
    % Find EEG samples that correspond to the beginning 
    % and end of the epochs of the current sub-task 
    sub_task_start_eeg = round(latency(...
        ismember(type,markers_sub_task_start)));
    sub_task_stop_eeg = round(latency(...
        ismember(type,markers_sub_task_stop)));
    
    % The first sub-task should begin where the 
    % task begin 
    if sub_task_order 
        sub_task_start_eeg(1) = task_start_eeg(1);
    end
    
    % No sub-task should start or stop outside
    % the limits defined for the task 
    sub_task_start_eeg(sub_task_start_eeg>...
        min(task_stop_eeg))=[];
    sub_task_stop_eeg(sub_task_stop_eeg>...
        min(task_stop_eeg))=[];

    % Define number of start and stop markers for the
    % current sub-task
    n_sub_task_start_eeg = length(sub_task_start_eeg);
    n_sub_task_stop_eeg = length(sub_task_stop_eeg);

    % Define a design matrix that is one for the samples that belong 
    % to the current sub-task, and zero for the samples that don't
    % This design matrix will be called sub_task_design
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