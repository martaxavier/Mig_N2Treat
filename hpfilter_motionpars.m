clear all
close all 

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','on');

fs=(1/1.26);
hpf=0.01;
 subjs = ["sub-32" "sub-35" "sub-36" "sub-37" "sub-38", ...
     "sub-39" "sub-40" "sub-42" "sub-43" "sub-44" "sub-45", ...
     "sub-46" "sub-47" "sub-48" "sub-49" "sub-50"]; 
task = 'task-rest';
dataset = 'NODDI';

for s = 1 : length(subjs) 
    
    % Specify input and output paths and data
    data_in = 'prefiltered_func_data_mcf.txt';
    data_out = 'prefiltered_func_data_mcf_tempfilt.txt';
    path_data_in = strcat(dataset,'/DATA/',subjs(s),'/',task,'/func/');
    path_data_out = path_data_in;
    path_img_out = strcat(path_data_in,'tempfilt');   
    if ~exist(path_img_out, 'dir'); mkdir(path_img_out); end

    % Perform HP filtering and save results
    rp=dlmread(fullfile(path_data_in,data_in));
    rp_hpf=highpass(rp,hpf,fs);
    dlmwrite(fullfile(path_data_out,data_out),rp_hpf);
    
    % Create vector of time-points
    n_pnts = size(rp,1);
    time = 0 : 1/fs : (n_pnts - 1) / fs;
    
    % Plot RP - Rotations - Original 
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Rotations)', 'FontSize',20);
    plot(time,rp(:,1:3)); xlabel('Time (sec)','FontSize',20); 
    ylabel('Rotations (rad)','FontSize',20); 
    legend('x','y','z'); img_out = 'RP_rotations';
    saveas(gcf,fullfile(path_img_out,img_out),'png');
    
    % Plot RP - Translations - Original 
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Translations)', 'FontSize',20);
    plot(time,rp(:,4:6)); xlabel('Time (sec)','FontSize',20);
    ylabel('Translations (mm)','FontSize',20); 
    legend('x','y','z'); img_out = 'RP_translations';
    saveas(gcf,fullfile(path_img_out,img_out),'png');
    
    % Plot RP - Rotations - HP filtered
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Rotations) - HP filtered',...
        'FontSize',20); plot(time,rp_hpf(:,1:3)); 
    xlabel('Time (sec)','FontSize',20); 
    ylabel('Rotations (rad)','FontSize',20);
    legend('x','y','z'); img_out = 'RP_rotations_tempfilt';
    saveas(gcf,fullfile(path_img_out,img_out),'png');
    
    % Plot RP - Translations - HP filtered  
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Translations) - HP filtered',...
        'FontSize',20); plot(time,rp_hpf(:,4:6)); 
    xlabel('Time (sec)','FontSize',20);
    ylabel('Translations (mm)','FontSize',20); 
    legend('x','y','z'); img_out = 'RP_translations_tempfilt';
    saveas(gcf,fullfile(path_img_out,img_out),'png');
    
end
 
