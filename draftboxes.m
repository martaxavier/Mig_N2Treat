subjects = ["sub-04" "sub-05" "sub-06" ...
	"sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12" ...
	"sub-13" "sub-14" "sub-15" "sub-16" "sub-17" "sub-19" ...
	"sub-20" "sub-21" "sub-22" "sub-23" "sub-24" "sub-25" ...
	"sub-27" "sub-29"];
sessions = ["run-1" "run-2" "run-3"];

n_subjects = length(subjects);
n_sessions = length(sessions);

nmse_test = zeros(n_subjects, n_sessions+1);
corr_test = zeros(n_subjects, n_sessions+1);
bic_test = zeros(n_subjects, n_sessions+1);

for s = 1 : n_subjects 

    for se = 1 : n_sessions
      
        data_in = strcat('lc4', '_', 'model', '_', 'blocked', '.mat');
		path_data_in = fullfile('PARIS\RESULTS', subjects(s), 'task-rest', sessions(se), 'models\ic_dmn_group\wavelet\l2_1');
		load(fullfile(path_data_in, data_in));
		nmse_test(s, se) = model.nmse;
        corr_test(s, se) = model.corr;
        bic_test(s, se) = model.bic;

    end

    data_in = strcat('lc4', '_', 'model', '_folds_', 'sessions', '.mat');
	path_data_in = fullfile('PARIS\RESULTS', subjects(s), 'task-rest\models\ic_dmn_group\wavelet\l2_1');
	load(fullfile(path_data_in, data_in));
    corr_test(s, 4) = mean(optimal.corr_test);
	nmse_test(s, 4) = mean(optimal.nmse_test);
    bic_test(s, 4) = mean(optimal.bic_test);
	
end

figure;
set(groot, 'defaultFigureUnits','centimeters');
set(groot, 'defaultFigurePosition', [0 0 22 30]);
b = boxplot(nmse_test, 'Colors', [0 0.4470 0.7410]);
set(b, {'linew'}, {2}); 
ax = gca; 
ax.XGrid = 'off'; ax.YGrid = 'on';
ax.YAxis.FontSize = 16;
set(gca,'XTickLabel', {'run-1', 'run-2', 'run-3', 'test-retest'}, 'FontSize', 22);
ylabel('NMSE Test', 'FontSize', 22);

figure;
set(groot, 'defaultFigureUnits','centimeters');
set(groot, 'defaultFigurePosition', [0 0 22 30]);
b = boxplot(corr_test, 'Colors', [0 0.4470 0.7410]);
set(b, {'linew'}, {2}); 
ax = gca; 
ax.XGrid = 'off'; ax.YGrid = 'on';
ax.YAxis.FontSize = 16;
set(gca,'XTickLabel', {'run-1', 'run-2', 'run-3', 'test-retest'}, 'FontSize', 22);
ylabel('Corr Test', 'FontSize', 22);

figure;
set(groot, 'defaultFigureUnits','centimeters');
set(groot, 'defaultFigurePosition', [0 0 22 30]);
b = boxplot(bic_test, 'Colors', [0 0.4470 0.7410]);
set(b, {'linew'}, {2}); 
ax = gca; 
ax.XGrid = 'off'; ax.YGrid = 'on';
ax.YAxis.FontSize = 16;
set(gca,'XTickLabel', {'run-1', 'run-2', 'run-3', 'test-retest'}, 'FontSize', 22);
ylabel('BIC Test', 'FontSize', 22);
