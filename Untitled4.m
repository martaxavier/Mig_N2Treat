% Mean
a = squeeze(nanmean(box_signal, 1));
a = squeeze(a(2,:));

% Max
b = squeeze(max(box_signal));
b = squeeze(b(2,:));
b = b - a;

% Min
c = squeeze(min(box_signal));
c = squeeze(c(2,:));
c = a - c;

% Error
d = max(b, c);

% Greater than 0.4
e = box_signal > 0.4;
e = squeeze(nansum(e));
e = squeeze(e(2,:));

% Total number
f = ~isnan(box_signal);
f = squeeze(nansum(f));
f = squeeze(f(2,:));

% Greater than 0.4 (%)
g = 100*e./f;

% Std
h = squeeze(nanstd(box_signal, 1));
h = squeeze(h(2,:));

%-----------------------------------------------

% 1-WAY ANOVAS
data_lc4_reliab = reshape(cell2mat(metric_reliab(1, 1)), [n_chans n_delays n_bands]);
data_rmsf_reliab = reshape(cell2mat(metric_reliab(2, 1)), [n_chans n_delays 1]);
data_tp_reliab = reshape(cell2mat(metric_reliab(3, 1)), [n_chans n_delays 1]);
data_wnd_reliab = reshape(cell2mat(metric_reliab(4, 1)), [n_chans n_delays n_bands]);

% CHANS 
data_reliab = cat(4, data_lc4_reliab, data_wnd_reliab);
g4 = zeros(size(data_reliab));
g4(logical(class_chans(:, 1)), :, :, :) = 1;
g4(logical(class_chans(:, 2)), :, :, :) = 2;
g4(logical(class_chans(:, 3)), :, :, :) = 3;
g4(logical(class_chans(:, 4)), :, :, :) = 4;
g4 = g4(:);

data_reliab_aux = cat(4, data_tp_reliab, data_rmsf_reliab);
g4_aux = zeros(size(data_reliab_aux));
g4_aux(logical(class_chans(:, 1)), :, :, :) = 1;
g4_aux(logical(class_chans(:, 2)), :, :, :) = 2;
g4_aux(logical(class_chans(:, 3)), :, :, :) = 3;
g4_aux(logical(class_chans(:, 4)), :, :, :) = 4;
g4_aux = g4_aux(:);
g4 = cat(1, g4, g4_aux);
data_reliab = cat(1, data_reliab(:), data_reliab_aux(:));

[p_c, tbl_c, signal_stats_c] = anovan(data_reliab(:), {g4(:)});
multcomp_c = multcompare(signal_stats_c); 


% DELAYS
data_reliab = cat(4, data_lc4_reliab, data_wnd_reliab);
g3 = zeros(size(data_reliab));
g3(:, 1, :, :) = 1;
g3(:, 2, :, :) = 2;
g3(:, 3, :, :) = 3;
g3(:, 4, :, :) = 4;
g3(:, 5, :, :) = 5;
g3(:, 6, :, :) = 6; 
g3 = g3(:);

data_reliab_aux = cat(4, data_tp_reliab, data_rmsf_reliab);
g3_aux = zeros(size(data_reliab_aux));
g3_aux(:, 1, :, :) = 1;
g3_aux(:, 2, :, :) = 2;
g3_aux(:, 3, :, :) = 3;
g3_aux(:, 4, :, :) = 4;
g3_aux(:, 5, :, :) = 5;
g3_aux(:, 6, :, :) = 6; 
g3_aux = g3_aux(:);
g3 = cat(1, g3, g3_aux);
data_reliab = cat(1, data_reliab(:), data_reliab_aux(:));

[p_d, tbl_d, signal_stats_d] = anovan(data_reliab(:), {g3(:)});
multcomp_d = multcompare(signal_stats_d); 


% BANDS
data_reliab = cat(4, data_lc4_reliab, data_wnd_reliab);
g2 = zeros(size(data_reliab));
g2(:, :, 1, :) = 1;
g2(:, :, 2, :) = 2;
g2(:, :, 3, :) = 3;
g2(:, :, 4, :) = 4;
g2  = g2(:);

data_reliab_aux = cat(4, data_tp_reliab, data_rmsf_reliab);
g2_aux = zeros(size(data_reliab_aux));
g2_aux(:, :, 1, :) = 5;
g2_aux(:, :, 1, :) = 5;
g2_aux(:, :, 1, :) = 5;
g2_aux(:, :, 1, :) = 5;
g2_aux = g2_aux(:);
g2 = cat(1, g2, g2_aux);
data_reliab = cat(1, data_reliab(:), data_reliab_aux(:));

[p_b, tbl_b, signal_stats_b] = anovan(data_reliab(:), {g2(:)});
multcomp_b = multcompare(signal_stats_b); 

% METRIC
data_reliab = cat(1, data_lc4_reliab(:), data_wnd_reliab(:), data_tp_reliab(:), data_rmsf_reliab(:));
g1 = zeros(size(data_reliab));

g1(1:1488) = 1;
g1(1489:2976) = 2;
g1(2977:3348) = 3;
g1(3349:3720) = 4;

[p_m, tbl_m, signal_stats_m] = anovan(data_reliab(:), {g1(:)});
multcomp_m = multcompare(signal_stats_m); 

% 4-way ANOVA

[p, tbl, signal_stats] = anovan(data_reliab(:), {g1(:), g2(:), g3(:), g4(:)}, 'model','linear');

% 4-way ANOVA
% data_lc4_reliab = reshape(cell2mat(metric_reliab(1, 1)), [n_chans n_delays n_bands]);
% data_rmsf_reliab = repmat(reshape(cell2mat(metric_reliab(2, 1)), [n_chans n_delays 1]), 1, 1, n_bands);
% data_tp_reliab = repmat(reshape(cell2mat(metric_reliab(3, 1)), [n_chans n_delays 1]), 1, 1, n_bands);
% data_wnd_reliab = reshape(cell2mat(metric_reliab(4, 1)), [n_chans n_delays n_bands]);
% 
% data_reliab = permute(cat(4, data_lc4_reliab, data_wnd_reliab, data_tp_reliab, data_rmsf_reliab), [4 1 2 3]);
% 
% g1 = zeros(size(data_reliab)); % metrics
% g2 = g1; % bands
% g3 = g1; % chans 
% g4 = g1; % delays
% 
% % lc4 and wnd
% for m = 1 : 2
%     g1(m, :, :, :) = m;
% end
% for m = 2 : n_metrics
%     g1(m, :, :, 1) = m;
% end
% for b = 1 : n_bands
%     g2(1:2, :, :, b) = b;
% end
% for c = 1 : size(class_chans, 2)
%     g3(1:2, logical(class_chans(:, c)), :, :) = c;
%     g3(3:4, logical(class_chans(:, c)), :, 1) = c;
% end
% for d = 1 : n_delays 
%     g4(1:2, :, d, :) = d;
%     g4(3:4, :, d, 1) = d;
% end
% 
%  g1(~g1) = NaN;
%  g2(~g2) = NaN;
% 
% [p, tbl, signal_stats] = anovan(data_reliab(:), {g1(:), g2(:), g3(:), g4(:)}, 'model','interaction');

% multcomp = multcompare(signal_stats); 