% Pre-allocate the matrix containing 
% the order of the acf model for each 
% subject 
order = zeros(length(subjects), 1);

% Only if some of the metrics specified
% are performed with the deconvolved BOLD
if contains(metrics, 'deconv')
    order_deconv = zeros(length(subjects), 1);
end

% ------------------------------------------------------------
% Go through subjects
% ------------------------------------------------------------
for s = 1 : length(subjects)
    
    % Create output directories, if non existent 
    if ~exist(char(path_img_out), 'dir'); mkdir(char(path_img_out)); end
    if ~exist(char(path_pars_out), 'dir'); mkdir(char(path_pars_out)); end
    
    % Load and read BOLD data of current subject 
    y = dlmread(char(fullfile(path_data_in(s), data_in)));
    
    % Estimate the order of the acf model 
    % for the BOLD time-series of the current 
    % subject 
    img_out = 'PAS_BOLD.png';
    order(s) = get_acf_order(y, acf_conf, path_img_out, img_out);
    
    if contains(metrics,'deconv')
            
        y_deconv = dlmread(char(fullfile(path_data_deconv_in(s), ...
            data_deconv_in)));
        img_out = 'PAS_BOLD_DECONV.png';
        order_deconv(s) = get_acf_order(y_deconv, acf_conf, ...
            path_img_out, img_out);
        
    end
    
end %finish looping through subjects 

% ------------------------------------------------------------
% Save output matrices with optimal parameters  
% ------------------------------------------------------------

% Create the output table 
subjs = subjects';
acf_order = table(subjs, order);

% Save table containing the estimated acf order for these
% datasets 
save(char(fullfile(path_pars_out, 'acf_order.mat')), 'acf_order');

% Repeat for deconvolved BOLD data 
if contains(metrics, 'deconv')
    acf_order = acf_order_deconv;
    acf_order_deconv = table(subjs, order_deconv);
    save(char(fullfile(path_pars_out, 'acf_order_deconv.mat')), ...
        'acf_order');
end


%/////////////////////////////////////////////////////////////
% SUBFUNCTIONS 
%/////////////////////////////////////////////////////////////

% ============================================================
% [order] = get_acf_order(y, conf, path_img_out, img_out)           
% ============================================================

function [order] = get_acf_order(y, conf, path_img_out, img_out)

%   [order] = get_arf_order(y) estimates the order 
%   of the autoregressive function of input signal y
%
%   INPUTS:
%
%     y                The signal y
%     conf             The confidence interval 
%   
%   OUTPUTS:
%
%     order            The order of the ARF

%------------------------------------------------
% Get order of the AR model
%------------------------------------------------

% Fit an AR(p) model using aryule
% to obtain the reflection coefficients 
p = 20; [~, ~, ref_coefs] = aryule(y, p);

% The negative of the reflection coefficients
% is the partial autocorrelation sequence
pacs = - ref_coefs;

% Estimate the confidence interval
n_pnts = length(y);
uconf = sqrt(2)*erfinv(conf)/sqrt(n_pnts);
lconf = -uconf;

% Check which is the first pacs value 
% to fall outside of the conf interval
outside = abs(pacs) < uconf;
outside = find(outside); 
order = outside(1);

%------------------------------------------------
% Plot the PACS and the confidence interval
%------------------------------------------------

stem(pacs, 'filled', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.7)
xlabel('Lag', 'FontSize', 16);
ylabel('Partial ACF', 'FontSize',16);
title('Partial Autocorrelation Sequence of Y',...
    'FontSize',16); xlim([1 p]); hold on
plot([1 p], [1 1]'*[lconf uconf], 'Color', ...
    [0.4660 0.6740 0.1880], 'LineWidth', 1);
grid on

saveas(gcf,char(fullfile(path_img_out, img_out)));




end