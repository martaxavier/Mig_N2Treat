function [conspec_sig, decision, p_values, p_thresh_corrected] = ...
    statistical_filtering(data, f_range, n_freq, f_vector, ...
    tf_sliding_win_seconds, fs_data, fs_cspec, metric, ...
    conspec, p_values, stat_filt_method, surr_method, ...
    n_surrs, path_img_out)

% Define some parameters 
% for the surrogate analysis 
win_cut_seconds = 10;
p_thresh_filt = 0.05;

% Get metric parameters
get_metric_pars
bands = bands;

% Read size of input data 
%n_chans = size(data,1);
n_pnts = size(data,2);

% Statistical filtering method 
switch stat_filt_method
    
    % Get p-values analytically 
    case 'analytical'
        
        if isempty(p_values)
            disp('Analytical p-values not available');
            decision = ones(size(conspec));
            p_thresh_corrected = 0;
            conspec_sig = conspec;
            return
            
        end
        
    % Get p-values through 
    % surrogate analysis 
    case 'surrogate'
        
        % Pre-allocate matrix of all connectivity surrogates 
        conspec = average_frequency(conspec, f_vector, bands, 3);
        conspec_surr = zeros([size(conspec) n_surrs]); 
        data_surr = zeros([size(data) n_surrs]);

        % Generate connectivity surrogates 
        parfor s = 1 : n_surrs

            data_par = data;
            bands_par = bands;

             % Generate the surrogate time-series 
             % All signals must be randomnize separately,
             % otherwise their connectivity remains 
             switch surr_method

                case 'phase_shuffle'

                    % Generate the current surrogates 
                    data_surr_par = phaseran(data_par', 1);
                    data_surr_par = data_surr_par';

                case 'block_shift'

                    % Create window in which the 
                    % time-series is restricted 
                    % to be cut 
                    mid = floor(n_pnts/2);
                    win_cut_samples = win_cut_seconds * fs_data + 1;
                    win_cut = mid - floor(win_cut_samples/2) : ...
                        mid + floor(win_cut_samples/2);

                    % Choose randomly one point from that window 
                    cut = round(win_cut(1) + rand(1)*...
                        (win_cut(end) - win_cut(1)));

                    % Cut the signal at that point and switch  
                    % the two resulting halves 
                    idxs_surr = [cut + 1 : n_pnts 1 : cut];
                    data_surr_par = data_par(:,idxs_surr);

             end   

            % Compute the connectivity spectrum 
            [conspec_surr_par, ~, ~] = tf_analysis_con_spectrum...
                (data_surr_par, f_range, n_freq, tf_sliding_win_seconds, ...
                fs_data, fs_cspec, metric, stat_filt_method);    
            data_surr(:, :, s) = data_surr_par;

            % Average connectivity for each frequency band 
            conspec_surr_par = average_frequency(conspec_surr_par, ...
                f_vector, bands_par, 3);   
            conspec_surr(:, :, :, s) = conspec_surr_par;

        end % finish looping through surrogates 

        % Compute the p-value - proportion of surrogates
        % that are higher than the original connectivity
        is_greater = conspec_surr > repmat(conspec, ...
            1, 1, 1, n_surrs);
        p_values = sum(is_greater, 4) ./ n_surrs;

        % No statistical filtering
        case 'none'

            p_values = zeros(size(conspec));
            decision = ones(size(conspec));
            p_thresh_corrected = 0;
            conspec_sig = conspec;
            return
            
end

% Correct for multiple comparisons, using the 
% FDR correction 
[decision, p_thresh_corrected, ~, ~] = ...
    fdr_bh(p_values, p_thresh_filt, 'pdep', 'no');
disp(size(p_values)); disp(p_thresh_corrected);

% Set to zero connectivity values that 
% didn't survive the significance test
conspec_sig = zeros(size(conspec));
conspec_sig(decision) = conspec(decision); 

% Leave if statistical filtering method 
% is not surrogate analysis 
if ~strcmp(stat_filt_method, 'surrogate')
    return
end

% Check that phase and amplitude dynamics 
% are preserved in the surrogates generated 
data_surr_X = fft(data);
data_X = fft(data);

% Check if the phase distribution is preserved
figure('Name', 'Distribution of data phase');
subplot(1, 2, 1); histogram(angle(data_X));
title('Original data'); xlabel('Phase (rad)');
subplot(1, 2, 2); histogram(angle(data_surr_X));  
title('Surrogate data'); xlabel('Phase (rad)');
img_out = strcat(upper(metric),'_phase_dist_surr.png'); 
saveas(gcf, fullfile(path_img_out, img_out));

% Check if the auto-correlation is preserved 
figure('Name', 'Auto-correlation of the data');
subplot(1, 2, 1); autocorr(double(data(1, :)));
title('Original data');
subplot(1, 2, 2); autocorr(double(data_surr(1, :, 1)));
title('Surrogate data');
img_out = strcat(upper(metric),'_auto_correlation_surr.png'); 
saveas(gcf, fullfile(path_img_out, img_out));

% Check if the mean distribution is preserved
figure('Name', 'Distribution of data mean');
subplot(1, 2, 1); histogram(mean(data, 2)); 
title('Original data'); xlabel('Mean');
subplot(1, 2, 2); histogram(mean(data_surr, 2)); 
title('Surrogate data'); xlabel('Mean');
img_out = strcat(upper(metric),'_mean_dist_surr.png'); 
saveas(gcf, fullfile(path_img_out, img_out));

% Check if the standard deviation is preserved 
figure('Name', 'Distribution of data standard deviation');
subplot(1, 2, 1); histogram(std(data, 0, 2)); 
title('Original data');  xlabel('Std');
subplot(1, 2, 2); histogram(std(data_surr, 0, 2)); 
title('Surrogate data');  xlabel('Std');
img_out = strcat(upper(metric),'_std_dist_surr.png'); 
saveas(gcf, fullfile(path_img_out, img_out));

end
        