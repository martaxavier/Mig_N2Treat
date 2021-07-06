% Performs temporal filtering and re-referencing of EEG data

%------------------------------------------------------------    
% Go through subjects
%------------------------------------------------------------ 

for s = 1 : length(subjects)
    
    % Broadcast the current pipeline stage 
    disp(strcat('Processing EEG for', ... 
    " ", subjects(s),' ...'));

    % Load input data matrix for current subject 
    EEG = pop_loadset(char(fullfile(path_data_in(s, se), data_in)));
    
    % Perform temporal filtering
    EEG_filtered = pop_eegfiltnew(EEG, ...
        highpass_filter, lowpass_filter);
    
    % Re-reference to average [] channel (exclude channel 32)
    EEG_avgref = pop_reref(EEG_filtered, [], 'exclude', 32);
    
    % Remove 32rd channel from data matrix
    EEG_avgref.data = EEG_avgref.data(1:31, :);
    EEG_avgref.chanlocs = EEG_avgref.chanlocs(1:31);
    EEG_avgref.nbchan = 31; EEG = EEG_avgref;
    
    % Save filtered and re-referenced EEG data
    if ~exist(path_data_out(s, se), 'dir'); mkdir(path_data_out(s, se)); end
    save(fullfile(path_data_out(s, se), data_out), 'EEG');

end

