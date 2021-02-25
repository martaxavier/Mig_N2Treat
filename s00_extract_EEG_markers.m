%------------------------------------------------------------    
% Go through subjects
%------------------------------------------------------------
        
for s = 1 : length(subjects)
    
    % Load eeg dataset of specified subject 
    load(fullfile(path_data_in(s),data_in));

    % Define current subject 
    subject = subjects(s);
    
    % Create output directories if non existent 
    if ~exist(path_data_eeg_out(s), 'dir')
        mkdir(path_data_eeg_out(s));
    end
    if ~exist(path_data_bold_out, 'dir')
        mkdir(path_data_bold_out); 
    end   
    
    disp(strcat('Extracting EEG markers', ...
        ' for subject'," ", subject, ' ...'));
 
    %---------------------------------------------------------    
    % Extract information from EEG event struct 
    %---------------------------------------------------------

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
    
    % For the NODDI dataset, remove 5 first scans 
    if strcmp(dataset, 'NODDI')
        idxs_scan(1 : 5) = [];
    end
    scan_eeg = round(latency(idxs_scan));
    first_eeg = round(latency(idxs_scan(1))); 
    last_eeg = round(latency(idxs_scan(end)));

    %---------------------------------------------------------    
    % Compute sub-task design 
    %---------------------------------------------------------

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
        if sub_task_order == 1
            sub_task_start_eeg(1) = task_start_eeg(1);
        end

        % No sub-task should start or stop outside
        % the limits defined for the task 
        sub_task_stop_eeg(sub_task_stop_eeg< ...
            sub_task_start_eeg(1))=[];
        sub_task_start_eeg(sub_task_start_eeg>...
            min(task_stop_eeg))=[];
        sub_task_stop_eeg(sub_task_stop_eeg>...
            min(task_stop_eeg))=[];

        % Now we must replace each of the sub-task start and
        % stop samples by the scan samples closest to them 
        % In case of a tie, the earlier scan sample is selected
        % (the same criteria is applied when cutting the BOLD 
        % time-series)
        dif = repmat(sub_task_start_eeg',[length(scan_eeg) 1])...
            - repmat(scan_eeg,[1 length(sub_task_start_eeg)]);
        [~,sub_task_start_bold] = min(abs(dif)); 
        sub_task_start_eeg = scan_eeg(sub_task_start_bold);

        dif = repmat(sub_task_stop_eeg',[length(scan_eeg) 1])...
            - repmat(scan_eeg,[1 length(sub_task_stop_eeg)]);
        [~,sub_task_stop_bold] = min(abs(dif)); 
        sub_task_stop_eeg = scan_eeg(sub_task_stop_bold);

        % The task should end in where the first task_stop
        % marker lies 
        [~,last_bold] = min(abs(task_stop_eeg(1)-scan_eeg));

        % Define number of start and stop markers for the
        % current sub-task
        n_sub_task_start_eeg = length(sub_task_start_eeg);
        n_sub_task_stop_eeg = length(sub_task_stop_eeg);

        % Define a 'sub_task_design' matrix that is one for the 
        % samples that belong to the current sub-task, and zero 
        % for the samples that don't
        design_table_eeg = cat(2,sub_task_start_eeg, ...
            ones(n_sub_task_start_eeg,1));
        design_table_eeg = cat(1,design_table_eeg, ...
            cat(2,sub_task_stop_eeg, zeros(n_sub_task_stop_eeg,1)));

        [~,order] = sort(design_table_eeg(:,1));
        design_table_eeg = design_table_eeg(order,:);
        n_markers = size(design_table_eeg,1);
        
        % Gain 1 TR per discontinuity (because upsampling/
        % downsampling is not additive)
        disc = find(diff(design_table_eeg(:,2))== - 1) + 1;
        design_table_eeg(disc) = design_table_eeg(disc) + ...
            mean(diff(scan_eeg));
        
        sub_task_design = zeros(1,n_pnts);

        for m = 1 : n_markers

            marker_start = design_table_eeg(m,1);

            if m == n_markers; marker_stop = n_pnts;
            else;marker_stop = design_table_eeg(m+1,1)-1;end

            marker_value = design_table_eeg(m,2);
            sub_task_design(marker_start:marker_stop) = marker_value;

        end

        % Define a 'design_table_bold' matrix that contains the  
        % sub-task design for the BOLD time-series 
        design_table_bold = cat(2,sub_task_start_bold, ...
            sub_task_stop_bold);
        design_table_bold = design_table_bold(order);  
        
        % THIS BIT OF CODE ONLY MAKES SENSE IF WE
        % THE LAST TASK IS ALSO THE FIRST TASK
        % (ACCOMODATE DIFFERENT POSSIBILITY)
        % The last sub-task should begin where the 
        % task ends 
        if sub_task_order == 1
            design_table_bold = cat...
                (2,design_table_bold,last_bold);
        else
            last_bold = sub_task_stop_bold(end);
        end     
        last_eeg = scan_eeg(last_bold);
        
        % Write the subject task timing file for the EEG time-series
        dlmwrite(fullfile(path_data_eeg_out(s), ...
            data_eeg_sub_task_out),sub_task_design);
        
        % Check if a timing file already exist and 
        % open one in case it doesn't
        if ~exist(fullfile(path_data_bold_out, ...
                data_bold_out), 'file')
        	fid = fopen(fullfile(path_data_bold_out, ...
                data_bold_out),'w');
        end
        
        % Don't generate file if already generated 
        if exist('fid','var')
            label = strcat(sub_task,'_',subject);
            fprintf(fid,'%s ',label);
            fprintf(fid,'%d ',design_table_bold); fprintf(fid,'\n');
        end
        
     end
     
     % Write the timing file for the EEG time-series
     first_last_eeg = [first_eeg last_eeg];
     dlmwrite(fullfile(path_data_eeg_out(s), ...
         data_eeg_out),first_last_eeg);
     
end

fclose('all');
