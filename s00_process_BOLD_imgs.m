% Performs upsampling and normalization of BOLD data
%

fs = fs_analysis;

%------------------------------------------------------------    
% Go through subjects
%------------------------------------------------------------ 

for s = 1 : length(subjects)
    
    disp(strcat('Processing BOLD images for', ...
        " ", subjects(s),'...'));
  
    BOLD = niftiread(fullfile(path_data_in(s),data_bold_in));
    DMN = niftiread(fullfile(path_data_in(s),data_dmn_in));
  
    %---------------------------------------------------------    
    % Binarize DMN mask and write dmn file 
    %--------------------------------------------------------- 
    
    DMN_bin = zeros(size(DMN));
    DMN_bin(DMN>=dmn_thr)=1;
    DMN_bin = reshape(DMN_bin, [numel(DMN) 1]);
    dlmwrite(fullfile(path_data_out(s),data_dmn_out),DMN_bin);
    
    %---------------------------------------------------------    
    % Upsample and normalize BOLD data
    %--------------------------------------------------------- 

    % Assign time vectors for upsampling 
    siz_BOLD = size(BOLD);
    n_pnts_bold = siz_BOLD(4); 
    time_bold = 0 : 1/fs_bold : (n_pnts_bold-1)/fs_bold;
    time = 0 : 1/fs : (n_pnts_bold-1)/fs_bold;
    n_pnts = length(time);
    
    % Save BOLD values that belong to the DMN 
    BOLD = reshape(BOLD, [prod(siz_BOLD(1:3)) ...
        n_pnts_bold]);
    BOLD = BOLD(DMN>=dmn_thr,:);
    
    % Pre-allocate object of interpolated values 
    BOLD_interp = zeros([size(BOLD,1) n_pnts]);
    BOLD_norm = BOLD_interp;
    
    parfor x = 1 : size(BOLD,1)
        
        % Upsample from fs_bold to fs 
        BOLD_interp(x,:) = interp1(time_bold, ...
            squeeze(BOLD(x,:)),time); 

        % Normalize BOLD data
        BOLD_norm(x,:) = zscore(BOLD_interp(x,:));

    end
    
    %---------------------------------------------------------    
    % Write BOLD output file 
    %--------------------------------------------------------- 
    
    dlmwrite(fullfile(path_data_out(s),data_bold_out),BOLD_norm');


end