% This script estimates the BOLD events and the HRF
% parameters of each subject's HRF and deconvolves
% their BOLD signal accordingly. 

n_subjects = length(subjects);

% ------------------------------------------------------------
% Write report in case of flag.report 
% ------------------------------------------------------------  

if flag.report == 2
    
    import mlflag.reportgen.dom.*;
    import mlflag.reportgen.flag.report.*;

    % Add section title 
    my_title = 'TIME-SERIES';
    H1 = get_report_heading(1,my_title);
    add(R,H1)
    
    for s = 1 : n_subjects
        
        % Add title for current subject 
        my_title = subjects(s);
        H2 = get_report_headinf(2,my_title);
        add(R,H2);
        
        % Add image of BOLD, BOLD deconvolved, BOLD
        % events
        source = fullfile(strcat(path_img_out(s, se), ...
            '\', rsn_method),strcat(data_bold_out, ...
            '.png'));
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
    end
   
    % Add section title 
    my_title = 'VOXEL-WISE';
    H1 = get_report_heading(1,my_title);
    add(R,H1)
    
    for s = 1 : n_subjects    
          
        % Add title for current subject 
        my_title = subjects(s);
        H2 = get_report_headinf(2,my_title);
        add(R,H2);
        
        % Load file with the results of deconvolution 
        deconv_file = strcat('deconv_',deconv_method,'.mat');
        load(strcat(path_data_out(s, se),deconv_file));
        
        % Compute the Pearson's correlation 
        % between heigth and the remaining 
        % parameters 
        val_height = hrf_pars(1,:);
        val_height(~val_height)=[]; 
        
        val_ttp = hrf_pars(2,:);
        val_ttp(~val_ttp)=[];
        [rho_ttp, pval_ttp] = ...
            corrcoef(abs(val_height),val_ttp);
        
        val_fwhm = hrf_pars(3,:);           
        val_fwhm(~val_fwhm)=[]; 
        [rho_fwmhm, pval_fwhm] = ...
            corrcoef(abs(val_height),val_fwhm);
        
        % Add histogram images
        my_title = 'Histograms';
        H3 = get_report_heading(3,my_title);
        add(R,H3);
        
        % Height
        source = fullfile(path_img_out(s, se),...
            'histo_estimated_hrfs_height.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Time-to-peak
        source = fullfile(path_img_out(s, se),...
            'histo_estimated_hrfs_timetopeak.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % FWHM
        source = fullfile(path_img_out(s, se),...
            'histo_estimated_hrfs_fwhm.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Number of events 
        source = fullfile(path_img_out(s, se),...
            'histo_estimated_number_events.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Add map images
        my_title = 'Maps';
        H3 = get_report_heading(3,my_title);
        add(R,H3);
        
        % Height
        source = fullfile(path_img_out(s, se),...
            'map_estimated_hrfs_height.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Time-to-peak
        source = fullfile(path_img_out(s, se),...
            'map_estimated_hrfs_timetopeak.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % FWHM
        source = fullfile(path_img_out(s, se),...
            'map_estimated_hrfs_fwhm.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % FWHM
        source = fullfile(path_img_out(s, se),...
            'map_estimated_number_events.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
    end

    return
    
end

% ------------------------------------------------------------
% Specify input parameters
% ------------------------------------------------------------  

% Input information 
fs = 4;
temporal_mask = [];

% Input basis function structure
xBF.name = 'Canonical HRF (time and dispersion derivatives)';  

xBF.len = 32;                            % length in seconds of basis
xBF.order = 1;                           % order of basis set

xBF.T = 16;                              % number of subdivisions of TR
                                         % microtime resolution 
xBF.TR = 1/fs;                           % scan repetition time (seconds)
xBF.T0 = fix(xBF.T/2);                   % first time bin 
xBF.dt = xBF.TR/xBF.T;                   % length of time bin (seconds)
                                         % xBF.dt = xBF.TR/xBF.T
                                         
xBF.AR_lag = 1;                          % AR(1) noise autocorrelation
xBF.thr = 1;                             % threshold for BOLD events 
 
min_onset_search = 2;                    % minimum delay allowed between
                                         % event and HRF onset (seconds)
max_onset_search = 10;                   % maximum delay allowed between 
                                         % event and HRF onset (seconds)

xBF.lag = fix(min_onset_search/xBF.dt)...% array of acceptable lags (bins)
    :fix(max_onset_search/xBF.dt);

% ------------------------------------------------------------
% Go through subjects 
% ------------------------------------------------------------  

best_delay = zeros(n_subjects,1);

for s = 1 : n_subjects
      
    if ~exist(path_data_out(s, se), 'dir')
        mkdir(path_data_out(s, se))
    end
    if ~exist(path_img_out(s, se), 'dir')
        mkdir(path_img_out(s, se))
    end

    data = dlmread(fullfile(path_data_in(s, se),data_in));
    n_pnts = size(data,1);
    n_voxs = size(data,2);
    
    if strcmp(deconv_method,'voxel_wise')
        dmn = dlmread(fullfile(path_data_in(s, se),dmn_in));
        data_img = zeros([size(dmn) 4]);
    end
    
    if flag.report == 0
    
        % ------------------------------------------------------------
        % Estimate BOLD events and HRF 
        % -------------------------------------------------------------  

        disp(strcat('Estimating HRF',...
        ' for', " ", subjects(s),' ...'));

        % Estimate the HRF basis function and BOLD events 
        % hrf_beta(1) - scale
        % hrf_beta(2) - onset
        % hrf_beta(3) - best lag 
        % bf - orthogonalized HRF basis function 
        % event_bold - estimated pseudo-events
        [hrf_beta, bf, event_bold] = ...
            rsHRF_estimation_temporal_basis(data,xBF,temporal_mask);

        if strcmp(deconv_method,'time_series')
            n_events=length(event_bold{1,1});
        end

        % Scale HRF basis function 
        % by the scaling parameter
        hrfa = bf*hrf_beta(1:size(bf,2),:);

        % Estimate HRF parameters from hrfa 
        % hrf_pars(1) - height
        % hrf_pars(2) - time to peak (derivative < 0.001)
        % hrf_pars(3) - width 
        hrf_pars = zeros(3,n_voxs);
        for v = 1 : n_voxs
            hrf_pars(:,v) = ...
                wgr_get_parameters(hrfa(:,v),xBF.TR/xBF.T);
        end 
    
    else
        
        % Load file with the results of deconvolution 
        deconv_file = strcat('deconv_',deconv_method,'.mat');
        load(fullfile(path_data_out(s, se),deconv_file));
        
    end
    
    num_events = zeros(length(event_bold),1);
    for i = 1 : length(event_bold)
        num_events(i) = length(event_bold{i});
    end
        
    if strcmp(deconv_method,'voxel_wise')
        data_img = squeeze(data_img);
        data_img(logical(dmn),1) = hrf_pars(1,:);
        data_img(logical(dmn),2) = hrf_pars(2,:);
        data_img(logical(dmn),3) = hrf_pars(3,:);
        data_img(logical(dmn),4) = num_events;
        data_img = reshape(data_img,...
            [n_xvoxs n_yvoxs n_zvoxs 4]);
    end

    % ------------------------------------------------------------
    % Retrieve estimated lag 
    % ------------------------------------------------------------  

    % Estimated HRF lag (seconds)
    %hrf_lag = hrf_beta(end)*xBF.dt;
    hrf_lag = hrf_pars(2,:); 

    % Save best lag for current subject 
    if strcmp(deconv_method,'time_series')
        best_delay(s) = hrf_lag;
    end

    % ------------------------------------------------------------
    % Perform deconvolution of the BOLD signal (default lag)
    % ------------------------------------------------------------  

    disp(strcat('Performing BOLD deconvolution',...
        ' for', " ", subjects(s),' ...'));

    hrfa_TR = resample(hrfa,1,xBF.T);
    hrf=hrfa_TR;

    % Deconvolution using a Wiener (restoration) filter 
    H=fft(cat(1,hrf, zeros(n_pnts-size(hrf,1),n_voxs)));
    B=fft(data);
    data_deconv = ifft(conj(H).*B./(H.*conj(H)+.1*mean(H.*conj(H))));

    % 0 mean and standard deviation 1 
    data_deconv = zscore(data_deconv);

    % Save deconvolved BOLD sinal in output directory
    dlmwrite(strcat(path_data_out(s, se),data_out,'.txt'), ...
        zscore(data_deconv));

    deconv_file = strcat('deconv_',deconv_method,'.mat');
    save(strcat(path_data_out(s, se),deconv_file),...
        'hrf_beta','hrf_pars','event_bold');
    
    % Leave if image generation
    % flag is off 
    if flag.report == 0 
        return
    end
    
    % ------------------------------------------------------------
    % Time-series: plot deconvolved data and save results
    % ------------------------------------------------------------  
    
    if strcmp(deconv_method,'time_series')
        
        % Plot HRF
        figure; 
        time = (1:length(hrfa(:,1)))*xBF.TR/xBF.T;
        plot(time,hrfa(:,1),'Color','#0072BD');
        delay_msg = strcat(upper('delay:'), " ", num2str(hrf_lag), ' s');
        text(time(end)-5,1,delay_msg,'FontSize', 12); hold on;
        xlabel('Time (s)'); title('Estimated HRF');
        saveas(gcf,fullfile(path_img_out(s, se),strcat('estimated_hrf.png')));

        event_plot=zeros(1,n_pnts);
        event_plot(event_bold{1,1})=1;
        event_plot(~event_plot)=nan;

        time = (1:n_pnts)*xBF.TR;

        % For the HRF with the default lag 
        figure;
        text(time(end)-50,-3.5,delay_msg,'FontSize', 12); hold on;
        plot(time,zscore(data),'LineWidth',0.7,'Color','#0072BD'); hold on;
        plot(time,zscore(data_deconv),'LineWidth',0.7,'Color','#77AC30');
        stem(time,event_plot,'Color','k','LineWidth',1);
        legend('BOLD','BOLD deconvolved','BOLD events','FontSize',14, ...
            'FontWeight','normal'); legend('boxoff');
        xlabel('Time (s)');
        saveas(gcf,fullfile(path_img_out(s, se),strcat(data_out,'.png')));
        
        % NOTE - zscore returns the z-score for each element of X 
        % such that columns of X are centered to have mean 0 and 
        % scaled to have standard deviation 1 
     
    % ------------------------------------------------------------
    % Voxel-wise: plot histograms and maps of HRF pars within the DMN 
    % ------------------------------------------------------------  
    
    elseif strcmp(deconv_method,'voxel_wise')
        
        % Load structural image of current subject 
        struct = niftiread(fullfile(path_struct_in(s),struct_in));
        struct = squeeze(struct(:,:,:,round(size(struct,4)/2)));
        
        % Create images for plotting 
        zslice = 35;
        img_height = squeeze(data_img(:,:,zslice,1));
        img_ttp = squeeze(data_img(:,:,zslice,2));
        img_fwhm = squeeze(data_img(:,:,zslice,3));
        img_be = squeeze(data_img(:,:,zslice,4));
        
        % Compute the Pearson's correlation between 
        % heigth and the remaining parameters 
        val_height = hrf_pars(1,:);
   
        val_ttp = hrf_pars(2,:);
        [rho_ttp, pval_ttp] = ...
            corrcoef(abs(val_height),val_ttp);
        
        val_fwhm = hrf_pars(3,:);           
        val_fwhm(~val_fwhm)=[]; 
        [rho_fwmhm, pval_fwhm] = ...
            corrcoef(abs(val_height),val_fwhm);
        
        % Plot all HRFs 
        fig = figure('Name','HRFs');
        fig.Position(3:4) = fig.Position(3:4)*5;
        xBF.bf = get_bf(xBF);
        plot(xBF.bf*hrf_beta(1:size(xBF.bf,2),:));
        title('Estimated HRFs of voxels in the DMN', ...
            'FontWeight','Normal'); set(gca,'FontSize',16); 
        colorbar; img_out = 'estimated_hrfs';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png');  
        
        % Plot histogram of height
        fig = figure('Name','Height - histogram');
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(hrf_pars(1,:));
        title('Height of the estimated HRFs of voxels in the DMN',...
        'FontWeight','Normal');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_hrfs_height';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png');  
             
        % Plot map of height 
        fig = figure('Name','Height - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        backgroundimage = struct(:,:,zslice);
        overlayimage = img_height;
        plot_map(backgroundimage, overlayimage);
        title(strcat('Height of the estimated HRFs of voxels', ...
            'in the DMN, axial slice (z ='," ", num2str(zslice),')'),...
            'FontWeight','Normal'); set(gca,'FontSize',16); colorbar; 
        img_out = 'map_estimated_hrfs_height';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png');  
        
        % Plot histogram of time-to-peak
        fig = figure('Name','Time-to-peak'); 
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(hrf_pars(2,:));
        title('Time-to-peak of the estimated HRFs of voxels in the DMN',...
        'FontWeight','Normal');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_hrfs_timetopeak';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png'); 
        
        % Plot map of time-to-peak 
        fig = figure('Name','Time-to-peak - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        backgroundimage = struct(:,:,zslice);
        overlayimage = img_ttp;
        plot_map(backgroundimage, overlayimage);
        title(strcat('Time-to-peak of the estimated HRFs of voxels', ...
            'in the DMN, axial slice (z ='," ", num2str(zslice),')'),...
            'FontWeight','Normal'); set(gca,'FontSize',16);
        img_out = 'map_estimated_hrfs_timetopeak';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png'); 
        
        % Plot weighted map of time-to-peak 
        % IM HERE, SCALE CONTINUES TO BE INCORRECT 
        %fig = figure('Name','Time-to-peak - map');
        %fig.Position(3:4) = fig.Position(3:4)*5;
        %backgroundimage = struct(:,:,zslice);
        %overlayimage = img_ttp.*zscore(abs(img_height));
        %plot_map(backgroundimage, overlayimage);
        
        % Plot histogram of FWHM
        fig = figure('Name','FWHM');
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(hrf_pars(3,:));
        title('FWHM of the estimated HRFs of voxels in the DMN',...
        'FontWeight','Normal');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_hrfs_fwhm';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png');  
        
        % Plot map of FWHM 
        fig = figure('Name','FWHM - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        backgroundimage = struct(:,:,zslice);
        overlayimage = img_fwhm;
        plot_map(backgroundimage, overlayimage);   
        title(strcat('FWHM of the estimated HRFs of voxels', ...
            'in the DMN, axial slice (z ='," ", num2str(zslice),')'),...
            'FontWeight','Normal'); set(gca,'FontSize',16);
        img_out = 'map_estimated_hrfs_fwhm';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png'); 

        % Plot histogram of number of events
        fig = figure('Name','# BOLD events - histogram');
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(num_events);
        title('# of BOLD events of voxels in the DMN',...
        'FontWeight','Normal');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_number_events';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png');  
        
        % Plot map of number of events  
        fig = figure('Name','# BOLD events - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        backgroundimage = struct(:,:,zslice);
        overlayimage = img_be;
        plot_map(backgroundimage, overlayimage);
        title(strcat('# of BOLD events of voxels', ...
            'in the DMN, axial slice (z ='," ", num2str(zslice),')'),...
            'FontWeight','Normal'); set(gca,'FontSize',16);
        img_out = 'map_estimated_number_events';
        saveas(gcf,fullfile(path_img_out(s, se),img_out),'png');      
        
    end

end % finish looping through subjects

if strcmp(deconv_method,'time_series')

    % Save table containing best_lag for each subject 
    subjects = subjects';
    deconv_delay = table(subjects,best_delay);
    save(strcat(path,'\PARS\',rsn_method,...
        '\deconv_delay.mat'),'deconv_delay');

end

% ============================================================
% [] = plot_map(backgroundimage, overlayimage)        
% ============================================================

function [] = plot_map(backgroundimage, overlayimage)

    % Make zeros NaNs in overlayimage 
    overlayimage(~overlayimage) = NaN;
    dimx = 1:100; dimy = dimx';

    % Plot the background image
    ax1 = axes;
    pcolor(ax1,dimx,dimy,backgroundimage); 
    set(gca, 'Color', 'black');
    colormap(ax1,'gray'); c1 = colorbar;
    shading flat;
    axis equal; hold on;

    % Plot the overlay image
    ax2 = axes;
    e1 = pcolor(ax2,dimx,dimy,overlayimage); 
    set(e1,'AlphaData',0);
    set(gca, 'Color', 'none'); 
    colormap(ax2,'parula'); colorbar;
    shading flat; axis equal;

    c1.Visible = 'off';
    ax1.Visible = 'off'; ax1.XTick = []; ax1.YTick = [];
    ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];      

end

% ============================================================
% [bf] = get_bf(xBF)     
% ============================================================

function [bf] = get_bf(xBF)

l = xBF.len;
T  =  xBF.T;
dt = xBF.dt;

[~,p]       = spm_hrf(dt,[],T);
p(end)      = l;
[bf,p]      = spm_hrf(dt,p,T);
    
% Add time derivative 
dp       = 1;
p(6)     = p(6) + dp;
D        = (bf(:,1) - spm_hrf(dt,p,T))/dp;
bf       = [bf D(:)];
p(6)     = p(6) - dp;

% Add dispersion derivative 
dp   = 0.01;
p(3) = p(3) + dp;
D    = (bf(:,1) - spm_hrf(dt,p,T))/dp;
bf   = [bf D(:)];

end

