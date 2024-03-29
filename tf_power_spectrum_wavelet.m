function [alpha, f_vector, wavelet_cf, power] = ...
    tf_power_spectrum_wavelet(data, f_range, n_freq, ...
    wav_seconds, fs, plot_TF)
%
%   function [alpha, vf, wavelet_cf, power] = tf_power_spectrum_wavelet 
%   (data, f_range, n_freq, siz_wavelet, fs, plot_TF) performs TF
%   decomposition of the input signal by Morlet wavelet convolution
%   in the time-domain (implemented by frequency-domain multiplication)
%
%   INPUTS:
%
%           data        - Input data (#time-points x #signals)
%           f_range     - Frequency range for the TF-decomposition 
%           n_freq      - Number of frequency bins for the TF-decomposition
%           wav_seconds - Number of seconds of Morlet wavelet  
%           fs          - Sampling frequency of the input data, 'X'
%           plot_TF     - 1 to plot, 0 to don't 
%
%   OUTPUTS:
%
%           alpha       -
%           f_vector    - Frequency vector used in thr tf-decomposition
%           wavelet_cf  - Convolution coefficients 
%           power       - Output power time-series of the input data  
%

% ------------------------------------------------------------
% Parameters of frequency decomposition 
% ------------------------------------------------------------ 

% Parameters for Wavelet convolution 
R = 7;                                      % Morlet Wavelet factor

% Build vector of frequencies 
f_max = f_range(2);                         % maximum frequency;
f_min = f_range(1);                         % minimum freuency;

alpha = (f_max/f_min)^(1/n_freq) - 1; 
f_vector = ((1+alpha).^(1:n_freq)) * f_min; % vector of frequencies 
                                            % (more frequency bins for
                                            % the lowest frquencies)
                                    
% ------------------------------------------------------------
% Get a wavelet for each frequency 
% ------------------------------------------------------------ 

% Compute the length of the convolution 
n_pnts = size(data, 1);
n_sigs = size(data, 2);
n_pnts_wavelet = wav_seconds*fs + 1;
n_conv = n_pnts + n_pnts_wavelet - 1;

% Pre-allocate matrix of Morlet wavelet coefficients
%wavelet_cf = zeros(n_freq, size(data));

n = 1;

% Pre-allocate matrix containing all wavelets
n_wav = wav_seconds * fs + 1;
wavelet = zeros(n_wav, n_freq);

% Go through all the frequencies 
for F = f_vector

    % Define Morlet wavelet parameters 
    UL = round(wav_seconds/2); LL = - UL; 
    std_f = F/R; std_t = 1/(2*pi*std_f); 
    FB = 2*(std_t^2); 
    
    % Compute the Morlet wavelet at frequency F
    [wavelet(:, n), ~] = cmorwavf(LL, UL, wav_seconds*fs + 1, FB, F);
    
    n = n + 1;

end

% ------------------------------------------------------------
% Perform time-domain convolution by freq.-domain multiplication
% ------------------------------------------------------------ 

% Compute the FFT along each column 
% The FFT will be padded with zeros to match the 
% length of the convolution 
waveletF = fft(wavelet, n_conv);
dataF = fft(data, n_conv);

waveletF = repmat(permute(waveletF, [1 3 2]), 1, n_sigs, 1);
dataF = repmat(dataF, 1, 1, n_freq);

% Perform time-domain convolution by frequency domain
% multiplication 
wavelet_cf = ifft(dataF .* waveletF, n_conv) ./ n_conv;

% Truncate the result of convolution to obtain a time-series
% that is the  same length as the original signal
wavelet_cf = wavelet_cf(floor(n_pnts_wavelet/2) ...
    : end - floor(n_pnts_wavelet/2) - 1, :, :);

% Compute power of convoluted signal
% Time though columns; frequency through lines 
power = abs(wavelet_cf) .^ 2; 

% % Plot power spectrum w/ TF + data
% if plot_TF
%     
%     % Colorscale
%     Dtf_f = log(abs(wavelet_cf) + .001); 
%     
%     figure('Name', 'TF')
%     
%     subplot(211)
%     imagesc((1 : size(power, 1)) ./ fs, ...
%         1: n_freq, squeeze(Dtf_f));
%     
%     Fplot = (log([1 2 5 8 13 20 30 45 50 70 90]) - ...
%         log(f_vector(1))) ./ log(1 + alpha) + 2;
%     hold on; plot([1 n_pnts],[Fplot', Fplot'],'k'); hold off
%     
%     set(gca,'YTick',Fplot)
%     set(gca,'YTickLabel', [1 2 5 8 13 20 30 45 50 70 90],'FontSize',12)
%     ylabel('Frequency (Hz)','FontSize', 26);
%     xlabel('Time (s)','FontSize',26);
%     colorbar;
%     
%     subplot(212)
%     plot((1 : size(power, 1)) ./ fs, data);
%     title('EEG Data','FontSize',16);
%     axis tight; grid minor;
%     ylabel('Amplitude (\muV)','FontSize', 26);
%     xlabel('Time (s)', 'FontSize', 26) 
%     
% end

end