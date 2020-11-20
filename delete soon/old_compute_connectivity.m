% Computes connectivity measures between all pairs of EEG electrodes,
% using the function FT_connectivityanalysis of the fieldtrip toolbox
%
% PIPELINE:
%
% 1. time-frequency decomposition for all N time-points (fs=250Hz);
%
% 2. remove NaN values from resulting fourier spectrum;
%
% 3. downsample the fourier spectrum to intermediate fs (fs=4Hz)
%
% 4. compute the connectivity measure for each time-point (fs=250Hz)
%    to obtain a CHANNEL (31) x CHANNEL(31) x BAND(5) x TIME(N) matrix 
%
% 5. average the fourier spectrum, in each time-point, across bands 
%    to obtain a CHANNEL(31) x BAND(5) x TIME(N) matrix 
% 
% 6. convolve all connectivity predictors with hrf (6 delays)
% 
% 7. cut all predictors from first to last sample 

close all
clear all

path='C:\Users\marta\Documents\LASEEB\MigN2Treat'; cd(path)
set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','on');

%---------------------------------------------------------------------%
%------------------------- Input information -------------------------%
%---------------------------------------------------------------------%

% Define anaysis type 
cfg.freqdec = 'No';                     % 'Yes', 'No'
cfg.downsamp = 'No';                    % 'Yes', 'No'
cfg.computecon = 'Yes';                  % 'Yes', 'No'
cfg.metric = 'Imcoh';                   % 'Imcoh','Wpli'
cfg.wnd = 'Yes';                        % 'Yes','No'
cfg.shift = 'Conv';                     % 'Conv', 'Delay'
cfg.plot = [3];                          % ids of the plots to be plotted 
cfg.data = 'eeg_preprocessed';  
cfg.subject = 'sub-patient006';  

% Define parameters for the
% frequency decomposition
n_freq = 100;                          % number of frequencies 
freq_range = [1 30];                   % frequency range 
R = 7;                                 % number of cycles      

% Define other parameters
fs_new = 4;                            % intermediate frequency (Hz)
delay_range = [10 8 6 5 4 2];          % predictor delasys (sec)
size_kern = 32;                        % size of hrf kernel (seconds)
n_delays = length(delay_range);

bands(1).name = 'Delta'; bands(1).range = [1 4];  % delta band range
bands(2).name = 'Theta'; bands(2).range = [4 8];  % theta band range
bands(3).name = 'Alpha'; bands(3).range = [8 13]; % alpha band range
bands(4).name = 'Beta'; bands(4).range = [13 30]; % beta band range
n_bands = size(bands,2);                                      

% Define plot variables 
plot_channel = [10 9];
plotting_delay = 3;             

% Specify directory where images should be saved 
% Define input and output data and images paths 
path_data_in=strcat(path,'\DATA\',cfg.subject,'\eeg\');
path_data_out=strcat(path,'\DERIVATIVES\',cfg.subject,'\eeg\');
path_img_out=strcat(path,'\RESULTS\',cfg.subject,'\eeg\connectivity');

%---------------------------------------------------------------------%
%------------------------- Load data (EEGLAB) ------------------------%
%---------------------------------------------------------------------%

% Load eeg dataset of specified subject 
load(strcat(path_data_in,cfg.data,'.mat'));

% Extract event information 
event = struct2cell(EEG.event);
type = event(1,:,:); type = squeeze(type);
latency = event(2,:,:); latency=squeeze(latency);
latency = cell2mat(latency);

% Map simultaneous fmri scan datapoints 
str=('Scan Start');idxScan = find(ismember(type,str));
idx_first = idxScan(1); first = round(latency(idx_first)); 
idx_second = idxScan(2); second = round(latency(idx_second));
idx_last =idxScan(end); last = round(latency(idx_last));

clear str idx_first idx_second idx_last

% Extract further data from dataset
data = EEG.data;                     % eeg data matrix 
n_chans = size(data,1);               % number of channels
n_pnts = size(data,2);                % number of datapoints
fs = EEG.srate;                      % sampling frequency (Hz)
TR = (second-first)/fs;              % bold repetition time (sec)
time_vector = 0:1/fs:(n_pnts-1)/fs;

%-------------------------------------------------------------------------%
%------------------- Parameters for frequency analysis -------------------%
%-------------------------------------------------------------------------%

% Build vector of frequencies
freq_min = freq_range(1); 
freq_max = freq_range(2);

alphaspctrm = (freq_max/freq_min)^(1/n_freq) - 1; 
freq_vector = ((1+alphaspctrm).^(1:n_freq)) * freq_min; % use more frequency
                                                        % bins for the lower
                                                        % frequencies
clear freq_min freq_max alphasctrm

% Assign frequency band indexes
for i = 1: n_bands
    bands(i).idx = find(freq_vector>=bands(i).range(1) ...
        & freq_vector<bands(i).range(2));
end


%-------------------------------------------------------------------------%
%--------------------- Compute connectivity metric -----------------------%
%---------------------- and downsample predictors ------------------------%

epoch_size = 2; % seconds
epoch_samples = epoch_size*fs +1;

time_vector_new = 0 : 1/fs_new : (n_pnts-1)/fs;
n_pnts_new = length(time_vector_new); 

cxx = zeros(n_chans,n_freq,n_pnts_new);
cxy = zeros(n_chans,n_chans,n_freq,n_pnts_new);
im_coh = zeros(n_chans,n_chans,n_freq,n_pnts_new);

offset = (epoch_samples-1)/2;
n_pnts_final = time_vector_new(end)*fs+1;
vector = round(1+offset:fs/fs_new:n_pnts_final-offset);

% Compute auto-cross-spectrum for each channel 
for c = 1 : n_chans
    
    x = data(c,:);
    t = 1+(n_pnts_new-length(vector))/2;

    for i = vector
        
        epoch = i-(epoch_samples-1)/2:i+(epoch_samples-1)/2;
        xepoch = x(epoch);
        cxx(c,:,t) = cpsd(xepoch,xepoch,hamming(250),1,freq_vector,fs);      
        t = t+1;
        
    end
    
end

cxx(:,:,1:(epoch_samples-1)/2)=1;

% Compute the cross-spectrum
for c1 = 1 : n_chans
    
    for c2 = 1 : c1
        
        if c2 == c1; continue; end
        
        t = 1+(n_pnts_new-length(vector))/2;

        for i = vector

            x = data(c1,:);y=data(c2,:);
            epoch = i-(epoch_samples-1)/2:i+(epoch_samples-1)/2;
            xepoch = x(epoch);yepoch=y(epoch);

            cxy(c1,c2,:,t) = cpsd(xepoch,yepoch,hamming(250),1,freq_vector,fs);  
            im_coh(c1,c2,:,t)=imag(squeeze(cxy(c1,c2,...
                :,t))'./squeeze(sqrt(cxx(c1,:,t).*cxx(c2,:,t))));
            
            t=t+1;

        end
    
    end
  
end


conspctrm = im_coh+permute(im_coh,[2,1,3,4]);

% Update first and last samples to the new frequency
first_new = round(first*fs_new/fs);
last_new = round(last*fs_new/fs); 
%ceil(first*fs_new/fs);
%floor(last*fs_new/fs); 

%-------------------------------------------------------------------------%
%-------------------- Average across frequency bands ---------------------%
%-------------------------------------------------------------------------%

% Average connectivity spectrum across frequency bands
conspctrm_avg = zeros(n_chans,n_chans,n_bands,n_pnts_new);

for i = 1 : n_bands
    conspctrm_avg(:,:,i,:) = mean(conspctrm(:,:,bands(i).idx,:),3); 
end 

%-------------------------------------------------------------------------%
%---------------------- Plot connectivity measures  ----------------------%
%-------------------------------------------------------------------------% 

% Save metric name for plotting
switch cfg.metric 
    case 'Imcoh'
        s = 'Imaginary Part of Coherency';
        ss = 'ipc';
    case 'Wpli'
        s = 'Weighted Phase Lag Index';
        ss = 'wpli';
end 

% Save channel label string for plotting 
chanlocs = struct2cell(EEG.chanlocs);
labels = string(squeeze(chanlocs(1,:,:)));

if ~isempty(find(cfg.plot==1))
    
    % Define baseline interval
    % baseline is between -baseline(1)
    % and -baseline(2), in seconds 
    baseline = [7.5 5]; 

    % Plot connectivity in all frequency bins between channel 
    % c3 (chan 5) and all other channels, throughout time

    figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
        ' and all other channels (part I), baseline subtracted');
    figure('Name',figure_name);

    n =1;
    for chan = 1 :round(n_chans/2)

        % Subtract baseline, consisting of imc time-averaged 
        % in the interval [-7.5,-5] seconds  
        idx = first_ne-(baseline(1)*fs_new):first_ne-(baseline(2)*fs_new);
        subtract = mean(squeeze(conspctrm(plotting_channel(1),chan,:,idx)),2);
        signal = squeeze(conspctrm(plotting_channel(1),chan,:,first_ne:last_new))-...
            subtract;

        subplot(4,4,n); imagesc(time_vector_new(first_ne:last_new),...
            freq_vector,signal); 
        colorbar;caxis([min(min(signal)) max(max(signal))]); 
        title(strcat(labels(plotting_channel(1)),",",labels(chan)));
        xlabel('Time (s)'); ylabel(ss); 

        n = n +1;

    end

    figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
        ' and all other channels (part II), baseline subtracted');
    figure('Name',figure_name);

    n = 1;
    for chan = round(n_chans/2)+1 : n_chans

        % Subtract baseline, consisting of imc time-averaged 
        % in the interval [-7.5,-5] seconds  
        idx = first_ne-(baseline(1)*fs_new):first_ne-(baseline(2)*fs_new);
        subtract = mean(squeeze(conspctrm(plotting_channel(1),chan,:,idx)),2);
        signal = squeeze(conspctrm(plotting_channel(1),chan,:,first_ne:last_new))-...
            subtract;

        subplot(4,4,n); imagesc(time_vector_new(first_ne:last_new),...
            freq_vector,signal); 
        colorbar; caxis([min(min(signal)) max(max(signal))]); 
        title(strcat(labels(plotting_channel(1)),",",labels(chan)));
        xlabel('Time (s)'); ylabel(ss); 

        n = n +1;

    end
    
    n_pnts_new = length(first_ne:last_new);
    timeaux = 0:1/fs_new:(n_pnts_new-1)/fs_new;
    idx = first_ne-(baseline(1)*fs_new):first_ne-(baseline(2)*fs_new);
    subtract = mean(squeeze(conspctrm(plotting_channel(1),plotting_channel(2),:,idx)),2);
    signal = squeeze(conspctrm(plotting_channel(1),plotting_channel(2),:,first_ne+150:last_new))-...
        subtract;
    figure;imagesc(timeaux(1:end-150),freq_vector(1:end-150),signal); 
    colorbar; caxis([min(min(signal)) max(max(signal))]); 
    %title(strcat('Imaginary Part of Coherency between channel', " ",...
     %   labels(plotting_channel(1))," and channel"," ",labels(plotting_channel(2))),...
      %  'FontSize',16);
    xlabel('Time (s)','FontSize',26); ylabel('Frequency (Hz)','FontSize',26); 
        
    for chan = 1 : n_chans

        % Plot topography of connectivty in every frequency band,
        % between all pairs of channels, averaged throughout time 
        figure_name = strcat(s, ' (time avg) for all pairs of channel',...
            " ", num2str(chan)); figure('Name',figure_name);

        for i = 1 : n_bands

            subplot(2,ceil(n_bands/2),i); 
            title(strcat(bands(i).name, ' band')); 
            topoplot(conspctrm_avg(chan,:,i,first_ne:last_new),EEG.chanlocs,...
                'electrodes','labels','whitebk','on', 'gridscale',100); 
            colorbar; caxis([min(min(conspctrm_avg(chan,:,i,first_ne:last_new)))...
               max(max(conspctrm_avg(chan,:,i,first_ne:last_new)))]);

        end

    end 


    % Plot connectivity in all frequency bands between channel 
    % c3 (chan 5) and all other channels, throughout time

    figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
        ' and all other channels (part I), baseline subtracted');
    figure('Name',figure_name);

    n =1;
    for chan = 1 :round(n_chans/2)

        % Subtract baseline, consisting of imc time-averaged 
        % in the interval [-7.5,-5] seconds  
        idx = first_ne-(baseline(1)*fs_new):first_ne-(baseline(2)*fs_new);
        subtract = mean(squeeze(conspctrm_avg(plotting_channel(1),chan,:,idx)),2);
        signal = squeeze(conspctrm_avg(plotting_channel(1),chan,:,first_ne:last_new))-...
            subtract;

        subplot(4,4,n); imagesc(time_vector_new(first_ne:last_new),...
            mean(range,2),signal); 
        colorbar; caxis([min(min(signal)) max(max(signal))]); 
        title(strcat(labels(plotting_channel(1)),",",labels(chan)));
        xlabel('Time (s)'); ylabel(ss); 

        n = n +1;

    end

    figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
        ' and all other channels (part II), baseline subtracted');
    figure('Name',figure_name);

    n = 1;
    for chan = round(n_chans/2)+1 : n_chans

        % Subtract baseline, consisting of imc time-averaged 
        % in the interval [-7.5,-5] seconds  
        idx = first_ne-(baseline(1)*fs_new):first_ne-(baseline(2)*fs_new);
        subtract = mean(squeeze(conspctrm_avg(plotting_channel(1),chan,:,idx)),2);
        signal = squeeze(conspctrm_avg(plotting_channel(1),chan,:,first_ne:last_new))-...
            subtract;

        subplot(4,4,n); imagesc(time_vector_new(first_ne:last_new),...
            mean(range,2),signal); 
        colorbar; caxis([min(min(signal)) max(max(signal))]); 
        title(strcat(labels(plotting_channel(1)),",",labels(chan)));
        xlabel('Time (s)'); ylabel(ss); 

        n = n +1;

    end
    
end
    
%---------------------------------------------------------------------%
%------------------- Compute Weighted Node Degree --------------------%
%---------------------------------------------------------------------%

if strcmp(cfg.wnd,'Yes')
    wnd_pred = squeeze(mean(conspctrm_avg,2));
    wnd_pred = permute(wnd_pred,[3,1,2]);
end

%---------------------------------------------------------------------%
%------------------- Remove redundant values and ---------------------%
%------------------ linearize connectivity matrix --------------------%

% Convert all 2D matrices (chan x chan) into upper 
% triangular matrices and save in 1D vectors the 
% linear indices of matrices' non-zero values 
aux = reshape(conspctrm_avg.*triu(ones(n_chans),1),...
    [n_chans*n_chans,n_bands,n_pnts_new]);
idx = find(aux(:,1,1)); aux = aux(idx,:,:); clear idx;
aux = permute(aux,[3,1,2]); pred.pred = aux; clear aux;

% Create pred field map, which maps each entry of the 1D
% vector to the column (pred.map(1)) and row (pred.map(2))
% indexes of the original (chan x chan) 2D matrix 
[columns,rows]=meshgrid(1:n_chans,1:n_chans); 
columns = (reshape(triu(columns,1),[n_chans*n_chans,1]));
rows = (reshape(triu(rows,1),[n_chans*n_chans,1]));
[~,~,columns]=find(columns); [~,~,rows]=find(rows);
pred.map = table(rows,columns); map=pred.map;

save('ConnectivityChansMap.mat','map');

map = table2array(map);


%---------------------------------------------------------------------%
%------------------- Mirror padding the signal -----------------------%
%----------------------- before convolution --------------------------%

% Mirror padd the signal at a length
% equal to kernal size + 1 (seconds)
aux = pred.pred(first_ne:last_new,:,:);

% Pre-direction 
padsize = max(first_ne-1,size_kern*fs_new);
aux = padarray(aux,padsize,...
    'symmetric','pre');

% Post-direction
padsize = max(n_pnts_new-last_new,size_kern*fs_new);
pred.pred = padarray(aux,padsize,...
    'symmetric','post');

%---------------------------------------------------------------------%
%----------------- Convolution with hrf or delay  --------------------%
%---------------------------------------------------------------------%

% Allocate predDelay matrix 
n_preds = size(pred.pred,2);
pred_delay = zeros(n_pnts_new,...
    n_preds,n_delays,n_bands);

if strcmp(cfg.wnd,'Yes')
    wnd_pred_delay = zeros(n_pnts_new,...
        n_chans,n_delays,n_bands);
end

switch cfg.shift 
    
    case 'Conv'
        
        plot_shift = 'Conv';

        pred_delay = ConvolveFeatures(pred.pred,...
            fs_new,delay_range,size_kern);
        
        pred_delay = pred_delay(first_ne:last_new,:,:,:);      
        
        % Cut the eeg predictor matrices to   
        % match the bold acquisition times 
        pred.pred=pred.pred(first_ne:last_new,:,:);
        
        % Weighted node degree
        if strcmp(cfg.wnd,'Yes')
            wnd_pred_delay = ConvolveFeatures(wnd_pred,...
            fs_new,delay_range,size_kern);       
            wnd_pred_delay = wnd_pred_delay(first_ne:last_new,:,:,:);
            wnd_pred = wnd_pred(first_ne:last_new,:,:);
        end
        
        % Update other variables 
        n_pnts_new = length(first_ne:last_new);
        time_vector_new = 0: 1/fs_new : (n_pnts_new-1)/fs_new;
            
    case 'Delay'
        
        plot_shift = 'delayed';
        pred_delay = DelayFeatures...
            (conspctrm_avg,fs_new,delay_range,[first_ne,last_new]);  
        
end        


%-------------------------------------------------------------------------%
%---------------- Normalize predictors (0 mean, 1 std) -------------------%
%-------------------------------------------------------------------------%

% Normalize predictor to zero mean and standard deviation one
pred_norm = pred.pred - (repmat(mean(pred.pred), n_pnts_new, 1)); 
pred_norm = pred_norm./(repmat(std(pred_norm), n_pnts_new, 1));

% Normalize delayed pred to 0 mean and standard deviation 1
pred_delay_norm = reshape(pred_delay, ...
    [n_pnts_new, n_preds*length(delay_range)*size(pred_delay,4), 1]);
pred_delay_norm = pred_delay_norm - ...
    (repmat(mean(pred_delay_norm), n_pnts_new, 1)); %demeaned
pred_delay_norm = pred_delay_norm./...
    (repmat(std(pred_delay_norm), n_pnts_new, 1)); %one std
pred_delay_norm = reshape(pred_delay_norm, ...
    [n_pnts_new, n_preds, length(delay_range),size(pred_delay,4)]);

if strcmp(cfg.wnd,'Yes')
    
    % Normalize predictor to zero mean and standard deviation one
    wnd_pred_norm = wnd_pred - (repmat(mean(wnd_pred), n_pnts_new, 1)); 
    wnd_pred_norm = wnd_pred_norm./(repmat(std(wnd_pred_norm), n_pnts_new, 1));

    % Normalize delayed pred to 0 mean and standard deviation 1
    wnd_pred_delay_norm = reshape(wnd_pred_delay, ...
        [n_pnts_new, n_chans*length(delay_range)*size(wnd_pred_delay,4), 1]);
    wnd_pred_delay_norm = wnd_pred_delay_norm - ...
        (repmat(mean(wnd_pred_delay_norm), n_pnts_new, 1)); %demeaned
    wnd_pred_delay_norm = wnd_pred_delay_norm./...
        (repmat(std(wnd_pred_delay_norm), n_pnts_new, 1)); %one std
    wnd_pred_delay_norm = reshape(wnd_pred_delay_norm, ...
        [n_pnts_new, n_chans, length(delay_range),size(wnd_pred_delay,4)]);
    
end


%-------------------------------------------------------------------------%
%---------------- Plot predictors (after normalization) ------------------%
%-------------------------------------------------------------------------%

if ~isempty(find(cfg.plot==2))
    
    delay = plotting_delay;
    map = table2array(pred.map);


    % Plot connectivity, for each frequency band, between  
    % channel c3 and all other channels, throughout time

    for i = 1 : n_bands

        figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
            ' and all other channels (part I), after HRF convolution',...
            ' (delay ',num2str(delay_range(delay)), ' sec) and normalization, (',...
            bands(i).name, " ", 'band)'); figure('Name',figure_name);

        n =1;

        for chan = 1 : round(n_chans/2)

            if chan == plotting_channel(1)
                continue
            end

            % Map plot channel pair from 2D matrix to 1D vector
            idx2D = sort([chan plotting_channel(1)]);
            idx1D = intersect(find(map(:,1)==idx2D(1)),...
                find(map(:,2)==idx2D(2)));

            signal = pred_delay_norm(:,idx1D,delay,i);

            subplot(4,4,n); plot(time_vector_new,signal); 
            title(strcat(labels(plotting_channel(1)),",",labels(chan)));
            xlabel('Time (s)'); ylabel('IPC'); 

            n = n + 1;

        end

        figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
            ' and all other channels (part II), after HRF convolution',...
            ' (delay ',num2str(delay_range(delay)), ' sec) and normalization, (',...
            bands(i).name, " ", 'band)'); figure('Name',figure_name);

        n =1;

        for chan = round(n_chans/2) + 1 : n_chans

            if chan == plotting_channel(1)
                continue
            end

            % Map plot channel pair from 2D matrix to 1D vector
            idx2D = sort([chan plotting_channel(1)]);
            idx1D = intersect(find(map(:,1)==idx2D(1)),...
                find(map(:,2)==idx2D(2)));

            signal = pred_delay_norm(:,idx1D,delay,i);

            subplot(4,4,n); plot(time_vector_new,signal); 
            title(strcat(labels(plotting_channel(1)),",",labels(chan)));
            xlabel('Time (s)'); ylabel(ss); 

            n = n + 1;

        end

    end

    % Plot connectivity, for all frequency bands and all 
    % delays, between channel c3 and fc1, throughout time

    % Map plot channel pair from 2D matrix to 1D vector
    idx2D = sort(plotting_channel);
    idx1D = intersect(find(map(:,1)==idx2D(1)),...
        find(map(:,2)==idx2D(2)));

    for i = 1 : n_bands 

    figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
        ' and channel', " ", labels(plotting_channel(2)),' after HRF convolution',...
        'and normalization, (', bands(i).name, " ", 'band)');
    figure('Name',figure_name);

        for d = 1:n_delays
            plot(time_vector_new,pred_delay_norm(:,idx1D,d,i)); hold on;
        end 

    title(figure_name); xlabel('Time (sec)'); ylabel(ss);
    legend(string(delay_range));


    end

end       

%-------------------------------------------------------------------------%
%---------------------------- Remove outliers ----------------------------%
%-------------------------------------------------------------------------%

% Threshold datapoints with value greater than mean+-3std
o = 3;

pred_norm(pred_norm>o)=o; pred_norm(pred_norm<-o)=-o;
pred_delay_norm(pred_delay_norm>o)=o; 
pred_delay_norm(pred_delay_norm<-o)=-o;

if strcmp(cfg.wnd,'Yes')
    wnd_pred_norm(wnd_pred_norm>o)=o; 
    wnd_pred_norm(wnd_pred_norm<-o)=-o;
    wnd_pred_delay_norm(wnd_pred_delay_norm>o)=o; 
    wnd_pred_delay_norm(wnd_pred_delay_norm<-o)=-o;
end 

%-------------------------------------------------------------------------%
%---------------- Plot predictors (after normalization -------------------%
%----------------------- and outlier removal) ----------------------------%

if ~isempty(find(cfg.plot==3))
    
    delay = plotting_delay;
        
    % Plot connectivity, for each frequency band, between  
    % channel c3 and all other channels, throughout time

        for i = 1 : n_bands

            figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
                ' and all other channels (part I), after HRF convolution',...
                ' (delay ',num2str(delay_range(delay)), ' sec) and normalization, (',...
                bands(i).name, " ", 'band)'); figure('Name',figure_name);

            n =1;

            for chan = 1 : round(n_chans/2)

                if chan == plotting_channel(1)
                    continue
                end

                % Map plot channel pair from 2D matrix to 1D vector
                idx2D = sort([chan plotting_channel(1)]);
                idx1D = intersect(find(map(:,1)==idx2D(1)),...
                    find(map(:,2)==idx2D(2)));

                signal = pred_delay_norm(:,idx1D,delay,i);

                subplot(4,4,n); plot(time_vector_new,signal); 
                title(strcat(labels(plotting_channel(1)),",",labels(chan)));
                xlabel('Time (s)'); ylabel('IPC'); 

                n = n + 1;

            end

            savename = strcat('IPC_',bands(i).name,num2str...
                (delay_range(plotting_delay)),...
                'sec',num2str(plotting_channel(1)),'chan_I.png');
            saveas(gcf,fullfile(path_img_out,savename));

            figure_name = strcat(s, ' between channel'," ",labels(plotting_channel(1)),...
                ' and all other channels (part II), after HRF convolution',...
                ' (delay ',num2str(delay_range(delay)), ' sec) and normalization, (',...
                bands(i).name, " ", 'band)'); figure('Name',figure_name);

            n =1;

            for chan = round(n_chans/2) + 1 : n_chans

                if chan == plotting_channel(1)
                    continue
                end

                % Map plot channel pair from 2D matrix to 1D vector
                idx2D = sort([chan plotting_channel(1)]);
                idx1D = intersect(find(map(:,1)==idx2D(1)),...
                    find(map(:,2)==idx2D(2)));

                signal = pred_delay_norm(:,idx1D,delay,i);

                subplot(4,4,n); plot(time_vector_new,signal); 
                title(strcat(labels(plotting_channel(1)),",",labels(chan)));
                xlabel('Time (s)'); ylabel(ss); 

                n = n + 1;

            end

            savename = strcat('IPC_',bands(i).name,num2str...
                (delay_range(plotting_delay)),...
                'sec',num2str(plotting_channel(1)),'chan_I.png');
            saveas(gcf,fullfile(path_img_out,savename));

        end
    
    % Plot weighted node degree for each channel and band 
    if strcmp(cfg.samples,'Entire')&& strcmp(cfg.wnd,'Yes')
            
         for i = 1 : n_bands 
        
             n = 1;
            
             for chan = 1 : n_chans
                 
             figure_name = strcat('Weighted node degree of channel'," ",...
                 labels(chan),' after HRF convolution (delay ',...
                 num2str(delay_range(delay)), ' sec) and normalization, (',...
                    bands(i).name, " ", 'band)'); figure('Name',figure_name);

                plot(time_vector_new, wnd_pred_delay_norm(:,chan,plotting_delay,i));
                title(figure_name); 
                xlabel('Time(s)'); ylabel('Amplitude');
                
                savename = strcat('WND_',bands(i).name,num2str...
                    (delay_range(plotting_delay)),...
                    'sec',num2str(chan),'chan.png');
                saveas(gcf,fullfile(path_img_out,savename));

            end   
            
         end 
    end 
    
end


%-------------------------------------------------------------------------%
%---------------------- Write eeg predictor files ------------------------%
%-------------------------------------------------------------------------%


% Build final connectivity predictors matrix
pred_norm = reshape(pred_norm, [n_pnts_new,...
    n_preds*size(pred_norm,3), 1]);
pred_delay_norm = reshape(pred_delay_norm, [n_pnts_new,...
    n_preds*n_delays*size(pred_delay_norm,4), 1]);

if strcmp(cfg.wnd,'Yes')
    wnd_pred_norm = reshape(wnd_pred_norm, [n_pnts_new,...
        n_chans*size(wnd_pred_norm,3), 1]);
    wnd_pred_delay_norm = reshape(wnd_pred_delay_norm, [n_pnts_new,...
        n_chans*n_delays*size(wnd_pred_delay_norm,4), 1]);
end

% Write txt file with 
% final predictors matrix      

dlmwrite(strcat('eeg_',ss,'_feature.txt'),...
    pred_norm);

dlmwrite(strcat('eeg_',ss,'_feature_',...
     cfg.shift,'.txt'),pred_delay_norm);

if strcmp(cfg.wnd,'Yes')
    dlmwrite(strcat('eeg_',ss,'_feature_wnd.txt'),...
        wnd_pred_norm);
    dlmwrite(strcat('eeg_',ss,'_feature',...
         '_wnd_',cfg.shift,'.txt'),wnd_pred_delay_norm);
end


