clear all
close all 

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','off');

fs=(1/1.26);
hpf=0.01;
 subjs = ["sub-02" "sub-04" "sub-05" "sub-06", ...
     "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12", ...
     "sub-13" "sub-14" "sub-15" "sub-16" "sub-17", "sub-19", ...
     "sub-20" "sub-21" "sub-22" "sub-23", "sub-24" "sub-25", ...
     "sub-27" "sub-29"]; 
task = 'task-rest';
run = 'run-3';
dataset = 'PARIS';

for s = 1 : length(subjs) 
    
    % Specify input and output paths and data
    data_in = 'prefiltered_func_data_mcf.txt';
    data_out = 'prefiltered_func_data_mcf_tempfilt.txt';
    path_data_in = fullfile(dataset, 'DATA', subjs(s), task, run, 'func');
    path_data_out = path_data_in;
    path_img_out = fullfile(path_data_in, 'tempfilt');   
    if ~exist(path_img_out, 'dir'); mkdir(path_img_out); end

    % Perform HP filtering and save results
    rp=dlmread(fullfile(path_data_in, data_in));
    rp_hpf=highpass(rp, hpf, fs);
    dlmwrite(fullfile(path_data_out, data_out), rp_hpf);
    
    % Create vector of time-points
    n_pnts = size(rp,1);
    time = 0 : 1/fs : (n_pnts - 1) / fs;
    
    % Plot RP - Rotations - Original 
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Rotations)', 'FontSize', 20);
    plot(time, rp(:, 1:3)); xlabel('Time (sec)', 'FontSize', 20); 
    ylabel('Rotations (rad)','FontSize', 20); 
    legend('x', 'y', 'z'); img_out = 'RP_rotations';
    saveas(gcf, fullfile(path_img_out, img_out), 'png');
    
    % Plot RP - Translations - Original 
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Translations)', 'FontSize', 20);
    plot(time,rp(:, 4:6)); xlabel('Time (sec)', 'FontSize', 20);
    ylabel('Translations (mm)', 'FontSize', 20); 
    legend('x', 'y', 'z'); img_out = 'RP_translations';
    saveas(gcf, fullfile(path_img_out, img_out), 'png');
    
    % Plot RP - Rotations - HP filtered
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Rotations) - HP filtered',...
        'FontSize', 20); plot(time, rp_hpf(:, 1:3)); 
    xlabel('Time (sec)', 'FontSize', 20); 
    ylabel('Rotations (rad)', 'FontSize', 20);
    legend('x', 'y', 'z'); img_out = 'RP_rotations_tempfilt';
    saveas(gcf,fullfile(path_img_out, img_out), 'png');
    
    % Plot RP - Translations - HP filtered  
    fig = figure(); fig.Position(3:4) = fig.Position(3:4)*5;
    title('Motion Realignment Parameters (Translations) - HP filtered',...
        'FontSize', 20); plot(time, rp_hpf(:, 4:6)); 
    xlabel('Time (sec)', 'FontSize', 20);
    ylabel('Translations (mm)', 'FontSize', 20); 
    legend('x', 'y', 'z'); img_out = 'RP_translations_tempfilt';
    saveas(gcf, fullfile(path_img_out, img_out), 'png');
    
end
 
