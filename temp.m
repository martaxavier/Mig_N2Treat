 subjects = ["sub-36" "sub-37" "sub-38" ...
     "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-45" ...
     "sub-46" "sub-47" "sub-48" "sub-49" "sub-50"];
path = 'C:\Users\marta\Documents\LASEEB\MigN2Treaty\NODDI\DATA';

for s = 1 : length(subjects)
    
    subject = subjects(s);
    path_in = fullfile(path,subject,'task-rest','eeg');
    load(fullfile(path_in,'eeg_processed.mat'));
    EEG = EEGICA;
    delete 'eeg_processed.mat';
    save('eeg_processed.mat','EEG');
    
end