function [ f_vector, power ] = ...
    tf_decomposition_welch(data, ic_or_channel, ...
    f_range, n_freq, fs, fs_new, plot_TF)
%
%   function [ vf, power, fs ] = tf_decomposition_welch(dataset,
%   data, ic_or_channel, filters, fs) performs time-frequency 
%   decomposition of the input EEG signal using the Welch method 
%
%   INPUTS:
%
%   data                 the EEG data matrix (channels x time)
%   ic_or_channel        the identification of the IC/channel in
%                        which the TF-decomposition will be done  
%   f_range              the frequency range for the TF-decomposition
%   fs                   the original sampling frequency of the data
%   plot_TF              1 to plot, 0 to don't 
%
%   OUTPUTS:
%
%   f_vector             the vector of frequency bins in which the 
%                        power will be evaluated 
%   power                the matrix that contains the EEG power-spectrum, 
%                        for each EEG time-point (time x frequency) 
%   fs_new               the output sampling frequency of the data 
%
%-------------------------------------------------------------------------%
% Prepare input for TF decomposition 
%-------------------------------------------------------------------------%

dataeeg = squeeze(data(ic_or_channel, :));

% Parameters for TF decomposition 
epoch_size = 2;                   % epoch size (seconds)
epoch_samples = epoch_size*fs + 1;% epoch size (samples)

% Build vector of frequencies 
fmax = f_range(2); % maximum frequency;
fmin = f_range(1); % minimum freuency;

alpha = (fmax/fmin)^(1/n_freq) - 1; 
f_vector = ((1+alpha).^(1:n_freq)) * fmin;  % vector of frequencies 
                                            % (more frequency bins for
                                            % the lowest frquencies)

%-------------------------------------------------------------------------%
% Perform TF decomposition using the Welch method 
%-------------------------------------------------------------------------%

% Compute the new vector of time-points 
n_pnts = length(dataeeg);
time_vector_new = 0 : 1/fs_new : (n_pnts-1)/fs;
n_pnts_new = length(time_vector_new); 

% Pre-allocate vector of power amplitude coefficients
cxx = zeros(n_freq,n_pnts_new);

offset = (epoch_samples-1)/2;
n_pnts_final = time_vector_new(end)*fs + 1;

% Vector containing the center of each epoch 
vector = round(1+offset : fs/fs_new : n_pnts_final-offset);

t = 1+(n_pnts_new-length(vector))/2;
x = dataeeg;
    
for i = vector

    epoch = i-(epoch_samples-1)/2:i+(epoch_samples-1)/2;
    xepoch = x(epoch);
    
    % Compute power spectral density for current epoch 
    % The cross-power spectral density (cpsd) of one signal  
    % with itself corresponds to its power spectral density 
    % Tapering is made with a Hamming window of 250 ms 
    % [] specifies 50% of overlap between adjacent segments 
    % f_vector specifies the frequencies in which the DFT 
    % is to be evaluated 
    cxx(:,t) = cpsd(xepoch,xepoch,hamming(250),[],f_vector,fs);
    
    t = t+1;

end

offset_size = round((1+offset)/(fs/fs_new));
power = cxx(:,1+offset_size:end-offset_size);
power = padarray(power, [0 offset_size], ...
    'replicate','both');

% Plot power spectrum w/ TF + data
if plot_TF
    
    % Colorscale
    Dtf_f = log(abs(power) + .001); 
    
    figure('Name', 'TF')
    
    %subplot(211)
    imagesc((1:size(power,2)) ./ fs, ...
        1:length(f_vector), squeeze(Dtf_f));
    
    Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) ...
        - log(f_vector(1))) ./ log(1 + alpha) + 2;
    hold on; plot([1 n_pnts_new],[Fplot', Fplot'],'k');
    hold off
    
    set(gca,'YTick',Fplot);
    set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90],'FontSize',12)
    ylabel('Frequency (Hz)','FontSize', 26);
    xlabel('Time (s)','FontSize',26);
    colorbar;
    
    %subplot(212)
    %plot((1:size(power,2)) ./ fs_new)
    %title('EEG Data','FontSize',16)
    %axis tight
    %grid minor
    %ylabel('Amplitude (\muV)','FontSize', 26);
    %xlabel('Time (s)', 'FontSize', 26) 
end


return;