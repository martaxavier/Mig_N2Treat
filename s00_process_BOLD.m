% Performs upsampling and normalization of BOLD data
%

%------------------------------------------------------------    
% Go through subjects
%------------------------------------------------------------ 

fs = fs_analysis;

for s = 1 : length(subjects)
    
    disp(strcat('Processing BOLD for', ... 
        " ", subjects(s),' ...'));
    
    % Create output directories, if non existent 
    if ~exist(path_img_out(s, se), 'dir'); mkdir(path_img_out(s, se)); end
    if ~exist(path_data_out(s, se), 'dir'); mkdir(path_data_out(s, se)); end
      
    % Load input data matrix for current subject 
    BOLD = load(fullfile(path_data_in(s, se), data_in));
    
    %---------------------------------------------------------    
    % Upsample BOLD data
    %--------------------------------------------------------- 

    % Assign time vectors for upsampling 
    n_pnts_bold = length(BOLD); 
    time_bold = 0 : 1/fs_bold : (n_pnts_bold-1)/fs_bold;
    time = 0 : 1/fs : (n_pnts_bold-1)/fs_bold;
    
    % Upsample from 1/TR to fs_new 
    BOLD = interp1(time_bold, BOLD, time); 
    %BOLD(end) = [];
    %time(end) = [];
    
    % Plot upsampled BOLD for current subject 
    figure('Name', strcat('Upsampled BOLD of ',...
        " ",subjects(s))); plot(time, BOLD);
    xlabel('Time(s)'); ylabel('Amplitude');
    title(strcat('BOLD data of', " ", subjects(s),...
        ' upsampled at', " ", num2str(fs), ' Hz'));

    %---------------------------------------------------------    
    % Normalize BOLD data 
    %---------------------------------------------------------
    
    % Normalize to zero mean and standard deviation one 
    BOLD_norm = zscore(BOLD); 

    % Plot normalized BOLD for current subject 
    figure('Name',strcat('Processed BOLD of',...
        " ",subjects(s))); plot(time,BOLD_norm);
    xlabel('Time(s)'); ylabel('Amplitude');
    title(strcat('BOLD data of', " ", subjects(s),...
        ' preprocessed',' (', num2str(fs), ' Hz, normalized)'));
    
    % Save current image in output image directory
    img_out = 'BOLD_preproc.png';
    saveas(gcf, fullfile(path_img_out(s, se), img_out));
    
    %---------------------------------------------------------    
    % Write BOLD output file  
    %---------------------------------------------------------
    
    dlmwrite(fullfile(path_data_out(s, se), data_out), BOLD_norm');

end