%------------------------------------------------------------    
% Go through subjects
%------------------------------------------------------------
        
for s = 1 : length(subjects)

    % Load eeg dataset of specified subject 
    load(fullfile(path_data_in(s, se), data_in));

    % Define current subject 
    subject = subjects(s);
    
    if exist('EEGCAP', 'var'); EEG = EEGCAP; end
    % Add chanlocs structure for current subject
    % Remove ecg label if still present 
    % Save result in directory containing the analysis parameters
    labels = struct2cell(EEG.chanlocs);
    ecg = find(strcmp(string(squeeze(labels(1, :, :))), 'ECG'));
    eog = find(strcmp(string(squeeze(labels(1, :, :))), 'EOG'));
    
    if ~isempty(ecg)
        labels(:, :, ecg) = [];
        EEG.data(ecg, :) = [];
        EEG.chanlocs = cell2struct(labels, fieldnames(EEG.chanlocs));
    end
    
    if ~isempty(eog)
        labels(:, :, eog) = [];
        EEG.data(eog, :) = [];
        EEG.chanlocs = cell2struct(labels, fieldnames(EEG.chanlocs));
    end
    
    % Initiate chanlocs structure
    % if it doesn't yet exist 
    if s == 1
        chanlocs = cell2struct(labels, fieldnames(EEG.chanlocs));
        sub_labels = strings(length(squeeze(labels(1, :, :))), ...
            length(subjects));
        sub_labels(:, s) = string(squeeze(labels(1, :, :)));
    else
        chanlocs(s, :) = cell2struct(labels, fieldnames(EEG.chanlocs));
        sub_labels(:, s) = string(squeeze(labels(1, :, :)));
    end
    
    save(fullfile(path_data_out(s, se), data_out), 'EEG');
    
    
end % subjects 

% Check if labels are the same for all subjects
% If so, save only one chanlocs structure containing 
% the common labels; if not, save a sub_chanlocs table
% containing the chanlocs structure for each subject 
% ADD POSSIBILITY OF RE-ORDERING LABELS 
same_labels = 1;
for l = 1 : size(sub_labels, 1)
    lab = sub_labels(l, 1);
    num_labels = find(strcmp(sub_labels, lab));
    if num_labels ~= length(subjects)
        same_labels = 0;
        break
    end
end

% Prepare and save outputs 
if same_labels
    
    chanlocs = chanlocs(1, :);
    save(fullfile(path_pars, 'chanlocs.mat', 'chanlocs'));
    
else
    
    subjs = subjects';
    sub_chanlocs = table(subjs, chanlocs);
    save(fullfile(path_pars, 'sub_chanlocs.mat'), 'sub_chanlocs');

end