fs = 250; 
data = EEG.data;
[n_chans, n_pnts] = size(data);
time = 0 : 1/fs : (n_pnts - 1)/fs;

% RAW EEG DATA 
figure; plot(time, data(:));
axis tight; grid minor;
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlabel('Time (s)', 'FontSize', 36) 
ylabel('Amplitude (\muV)','FontSize', 36);

% WAVELET 
time_wav = 0 : 1/fs : wav_seconds;
figure; plot(time_wav, real(wavelet), ...
    'LineWidth',2);
axis tight; grid minor;
ax = gca;
ax.XAxis.FontSize = 16;
xlabel('Time (s)', 'FontSize', 36) 

% TF DATA 
Dtf_f = log(abs(wavelet_cf) + .001); 
figure('Name', 'TF')
imagesc((1:size(power,2)) ./ fs, ...
f_vector, squeeze(Dtf_f));
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ylabel('Frequency (Hz)','FontSize', 36);
xlabel('Time (s)','FontSize', 36);
colorbar;

% TF DATA BANDS
feature_d = eeg_features(:,10,1); 
feature_t = eeg_features(:,10,2); 
feature_a = eeg_features(:,10,3); 
feature_b = eeg_features(:,10,4);

figure('Name','TF Bands');
% plot(time, zscore(feature_b), 'LineWidth', 0.6); hold on;
% plot(time, zscore(feature_a) + max(zscore(feature_b)), ...
%     'Color', '#4DBEEE', 'LineWidth', 0.6);
plot(time, zscore(feature_t), ...
    'Color', '#77AC30', 'LineWidth', 0.8); hold on;
plot(time, zscore(feature_d) + max(zscore(feature_t)), ...
    'Color', '#EDB120', 'LineWidth', 0.8);
axis tight; grid minor
ax = gca;
ax.XAxis.FontSize = 16;
xlabel('Time (s)','FontSize', 36);
set(gca,'ytick',[]);

% HRF
% Assign time vector with fs 250 Hz 
time_hrf = 0 : 1/fs : kern_seconds;

% Plot hrfs 

for d = 1 : size(hrf,2)
    figure;
    plot(time_hrf, hrf(:,d), 'LineWidth', 3); hold on 
    axis tight; grid minor
    ax = gca;
    ax.XAxis.FontSize = 20;
    xlabel('Time (s)','FontSize', 40);
    set(gca,'ytick',[]);
end

% legend("10 seconds", "8 seconds", "6 seconds",...
%     "5 seconds", "4 seconds", "2 seconds",'FontSize',14)

% TF DATA BANDS CONVOLVED 
feature_d = features_delayed(:,10,1,1); 
feature_t = features_delayed(:,10,2,1); 
feature_a = features_delayed(:,10,3,1); 
feature_b = features_delayed(:,10,4,1);

n_pnts = size(features_delayed,1);
fs = 250;
time = 0 : 1/fs : (n_pnts - 1)/fs;

figure('Name','TF Bands');
plot(time, zscore(feature_b), 'LineWidth', 1); hold on;
plot(time, zscore(feature_a) + max(zscore(feature_b)), ...
    'Color', '#4DBEEE', 'LineWidth', 1);
% plot(time, zscore(feature_t), ...
%     'Color', '#77AC30', 'LineWidth', 1.2); hold on;
% plot(time, zscore(feature_d) + max(zscore(feature_t)), ...
%     'Color', '#EDB120', 'LineWidth', 1.2);
axis tight; grid minor
ax = gca;
ax.XAxis.FontSize = 16;
xlabel('Time (s)','FontSize', 36);
set(gca,'ytick',[]);

% CONNECTIVITY 
signal = squeeze(cross_spectrum(:, :, 10, 9));
Dtf_f = log(abs(signal) + .001); 
figure('Name','Cross-spectrum');
imagesc((1:size(signal,2)) ./ fs, ...
    f_vector, squeeze(Dtf_f));
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ylabel('Frequency (Hz)','FontSize', 36);
xlabel('Time (s)','FontSize', 36);
colorbar;