% Assigns the parameters that define the current metric 
% To add a new model, add a new section the end of the  
% script

%-------------------------------------------------------------
% Load chanlocs structure 
%-------------------------------------------------------------

if exist('subject', 'var')
    
    if exist(fullfile(path_pars, 'sub_chanlocs.mat'), 'file')    
        load(fullfile(path_pars, 'sub_chanlocs.mat'));
        idx = find(ismember(table2array(sub_chanlocs(:, 1)), subject));
        chanlocs = table2array(sub_chanlocs(idx, 2));  
        labels = struct2cell(chanlocs);
        labels = labels(1, :, :); 
        labels = string(squeeze(labels));  
    else  
        chanlocs = table2array(sub_chanlocs(1, 2));  
        labels = struct2cell(chanlocs);
        labels = labels(1, :, :); 
        labels = string(squeeze(labels)); 
        labels = repmat("Undefined", size(labels));   
    end
        
elseif exist(fullfile(path_pars, 'chanlocs.mat'), 'file')
    
    load(fullfile(path_pars, 'chanlocs.mat'));
	labels = struct2cell(chanlocs);
	labels = labels(1, :, :); 
	labels = string(squeeze(labels));
    
end
    
%--------------------------------------------------------
% Write parameters for each model 
%--------------------------------------------------------

switch metric 
    
    case 'lc4'
        
        % ================================================
        % LC4 (lin comb of band-specific power, 4 bands)           
        % ================================================

        id =            'lc4';
        eeg_metric =    'lc4';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
    
    case 'lc4_delay'
        
        % ================================================
        % LC4_DELAY (lin comb of band-specific power, 4 bands)           
        % ================================================

        id =            'lc4_delay';
        eeg_metric =    'lc4';
        eeg_shift =     'delay';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
        
    case 'lc4_deconv'
        
        % ================================================
        % LC4_DECONV (lin comb of band-specific power, 4 bands)           
        % ================================================

        id =            'lc4_deconv';
        eeg_metric =    'lc4';
        eeg_shift =     '';
        bold_shift =    'deconv';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     [];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
        
    case 'lc5'
        
        % ================================================
        % LC5 (lin comb of band-specific power, 5 bands)           
        % ================================================

        id =            'lc5';
        eeg_metric =    'lc5';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30; 30 60]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta", "Gamma"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
        
    case 'lc5_deconv'
        
        % ================================================
        % LC5_DECONV (lin comb of band-specific power, 5 bands)           
        % ================================================

        id =            'lc5_deconv';
        eeg_metric =    'lc5';
        eeg_shift =     '';
        bold_shift =    'deconv';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30; 30 60]';
        delays =        [];
        id_chans =      labels;
        id_delays =     [];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta", "Gamma"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];        

    case 'lc6'
        
        % ================================================
        % LC6 (lin comb of band-specific power, 6 bands)           
        % ================================================

        id =            'lc6';
        eeg_metric =    'lc6';
        eeg_shift =     'conv';
        bold_shift =          '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 10; ...
                        10 12; 12 20; 20 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha1", ...
                        "Alpha2", "Beta1", "Beta2"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
      
    case 'lc6_delay'
          
        % ================================================
        % LC6_DELAY (lin comb of band-specific power, 6 bands)           
        % ================================================

        id =            'lc6_delay';
        eeg_metric =    'lc6';
        eeg_shift =     'delay';
        bold_shift =          '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 10; ...
                        10 12; 12 20; 20 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha1", ...
                        "Alpha2", "Beta1", "Beta2"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];        
        
    case 'lc6_deconv'
        
        % ================================================
        % LC6_DECONV (lin comb of band-specific power, 6 bands)           
        % ================================================

        id =            'lc6_deconv';
        eeg_metric =    'lc6';
        eeg_shift =     '';
        bold_shift =          'deconv';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 10; ...
                        10 12; 12 20; 20 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     [];
        id_bands =      ["Delta", "Theta", "Alpha1", ...
                        "Alpha2", "Beta1", "Beta2"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];     
        
    case 'tp'
        
        % ================================================
        % TP (total power across all frequencies)           
        % ================================================

        id =            'tp';
        eeg_metric =    'tp';
        eeg_shift =     'conv';
        bold_shift =          '';
        chans =         1:length(labels);
        bands =         [1 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["TP"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
  
        case 'tp_delay'
        
        % ================================================
        % TP_DELAY (total power across all frequencies)           
        % ================================================

        id =            'tp_delay';
        eeg_metric =    'tp';
        eeg_shift =     'delay';
        bold_shift =          '';
        chans =         1:length(labels);
        bands =         [1 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["TP"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];

    case 'tp_deconv'
        
        % ================================================
        % TP_DECONV (total power across all frequencies)            
        % ================================================

        id =            'tp_deconv';
        eeg_metric =    'tp';
        eeg_shift =     '';
        bold_shift =    'deconv';
        chans =         1:length(labels);
        bands =         [1 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     [];
        id_bands =      ["TP"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands]; 

    case 'rmsf'
        
        % ================================================
        % RMSF (root mean squared frequency)           
        % ================================================

        id =            'rmsf';
        eeg_metric =    'rmsf';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["RMSF"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];

    case 'rmsf_delay'
        
        % ================================================
        % RMSF_DELAY (root mean squared frequency)           
        % ================================================

        id =            'rmsf_delay';
        eeg_metric =    'rmsf';
        eeg_shift =     'delay';
        bold_shift =          '';
        chans =         1:length(labels);
        bands =         [1 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["RMSF"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
        
    case 'rmsf_deconv'
        
        % ================================================
        % RMSF_DECONV (root mean squared frequency)             
        % ================================================

        id =            'rmsf_deconv';
        eeg_metric =    'rmsf';
        eeg_shift =     '';
        bold_shift =    'deconv';
        chans =         1:length(labels);
        bands =         [1 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     [];
        id_bands =      ["RMSF"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands]; 
        
    case 'icoh'
        
        % ================================================
        % ICOH (Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh';
        eeg_metric =    'icoh';
        con_metric =    'icoh';
        net_metric =    [];
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands, 2);
        dim =           [n_chans, n_chans, n_delays, n_bands];      
        
    case 'icoh_deconv'
        
        % ================================================
        % ICOH DECONV. (Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh_deconv';
        eeg_metric =    'icoh';
        con_metric =    'icoh';
        net_metric =    [];
        eeg_shift =     'deconv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands, 2);
        dim =           [n_chans, n_chans, n_delays, n_bands];            

    case 'icoh_wnd'
        
        % ================================================
        % ICOH WND (WND of Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh_wnd';
        eeg_metric =    'icoh_wnd';
        con_metric =    'icoh';
        net_metric =    'wnd';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];  
        
    case 'icoh_wnd_deconv'
        
        % ================================================
        % ICOH WND DECONV. (WND of Imaginary Part of Coh.)           
        % ================================================

        id =            'icoh_wnd_deconv';
        eeg_metric =    'icoh_wnd';
        con_metric =    'icoh';
        net_metric =    'wnd';
        eeg_shift =     'deconv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];         
        
    case 'icoh_wne'
        
        % ================================================
        % ICOH WNE (WNE of Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh_wne';
        eeg_metric =    'icoh_wne';
        con_metric =    'icoh';
        net_metric =    'wne';        
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
 
    case 'icoh_bc'
        
        % ================================================
        % ICOH BC (BC of Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh_bc';
        eeg_metric =    'icoh_bc';
        con_metric =    'icoh';
        net_metric =    'bc';        
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
        
    case 'icoh_bc_deconv'
        
        % ================================================
        % ICOH BC DECONV. (WND of Imaginary Part of Coh.)           
        % ================================================

        id =            'icoh_bc_deconv';
        eeg_metric =    'icoh_bc';
        con_metric =    'icoh';
        net_metric =    'wnd';
        eeg_shift =     'deconv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      1;
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];           
        
    case 'wpli'
        
        % ================================================
        % WPLI (Weighted Phase Lag Index)           
        % ================================================

        id =            'wpli';
        eeg_metric =    'wpli';
        con_metric =    'wpli';
        net_metric =    [];
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_chans, n_delays, n_bands];        

    case 'wpli_wnd'
        
        % ================================================
        % WPLI WND (WND of Weighted Phase Lag Index)           
        % ================================================

        id =            'wpli_wnd';
        eeg_metric =    'wpli_wnd';
        con_metric =    'wpli';
        net_metric =    'wnd';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:length(labels);
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];          
        
end

%-------------------------------------------------------------
% Classify channel locations 
%-------------------------------------------------------------

id_class_chans = {'Frontal', 'Parietal', 'Temporal', 'Occipital'};
%id_class_chans = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'Midline'};
class_chans = zeros(length(id_chans), length(id_class_chans));

% class_chans(:, 1) = contains(id_chans, 'F');
% class_chans(:, 3) =  contains(id_chans, 'T');        
% class_chans(:, 4) = contains(id_chans, 'O');        
% class_chans(:, 2) = contains(id_chans, 'C');    
class_chans(:, 1) = contains(id_chans, 'F') ...
                    & ~contains(id_chans, 'C') ...
                    & ~contains(id_chans, 'T');    
class_chans(:, 3) =  contains(id_chans, 'T');        
class_chans(:, 4) = contains(id_chans, 'O');      
class_chans(:, 2) = ~sum(class_chans, 2);       
%class_chans(:, 5) = contains(id_chans, 'Z');

