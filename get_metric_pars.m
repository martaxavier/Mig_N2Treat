% Assigns the parameters that define the current metric 
% To add a new model, add a new section the end of the  
% script

%-------------------------------------------------------------
% Load chanlocs structure 
%-------------------------------------------------------------

load(fullfile('PARS','chanlocs.mat'));

labels = struct2cell(chanlocs);
labels = labels(1,:,:); 
labels = string(squeeze(labels));

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
        chans =         1:31;
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
        chans =         1:31;
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
        chans =         1:31;
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [];
        id_chans =      labels;
        id_delays =     [];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
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
        chans =         1:31;
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
        chans =         1:31;
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
        chans =         1:31;
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
        chans =         1:31;
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
        chans =         1:31;
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
        bold_shift =          'deconv';
        chans =         1:31;
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
        chans =         1:31;
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
        chans =         1:31;
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
        chans =         1:31;
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
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:31;
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_chans, n_delays, n_bands];        

    case 'icoh_wnd'
        
        % ================================================
        % ICOH WND (WND of Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh_wnd';
        eeg_metric =    'icoh_wnd';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:31;
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];  
        
    case 'icoh_wne'
        
        % ================================================
        % ICOH WNE (WNE of Imaginary Part of Coherency)           
        % ================================================

        id =            'icoh_wne';
        eeg_metric =    'icoh_wne';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:31;
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
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:31;
        bands =         [1 4; 4 8; 8 12; 12 30]';
        delays =        [2 4 5 6 8 10];
        id_chans =      labels;
        id_delays =     ["2", "4", "5", "6", "8", "10"];
        id_bands =      ["Delta", "Theta", "Alpha", "Beta"];
        n_chans =       length(chans);
        n_delays =      length(delays);
        n_bands =       size(bands,2);
        dim =           [n_chans, n_delays, n_bands];
        
    case 'wpli'
        
        % ================================================
        % WPLI (Weighted Phase Lag Index)           
        % ================================================

        id =            'wpli';
        eeg_metric =    'wpli';
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:31;
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
        eeg_shift =     'conv';
        bold_shift =    '';
        chans =         1:31;
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
