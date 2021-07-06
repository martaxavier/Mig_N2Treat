% Check if the mean distribution is preserved
idxs_delta = 3:5; idxs_theta = 6:9;
idxs_alpha = 10:13; idxs_beta = 14:31;

figure('Name', 'Distribution of ICoh values');
subplot(2, 2, 1); histogram(conspec(:,:,idxs_delta));
title('Delta'); xlabel('Connectivity (Icoh)');
subplot(2, 2, 2); histogram(conspec(:,:,idxs_theta)); 
title('Theta'); xlabel('Connectivity (Icoh)');
subplot(2, 2, 3); histogram(conspec(:,:,idxs_alpha));
title('Alpha'); xlabel('Connectivity (Icoh)');
subplot(2, 2, 4); histogram(conspec(:,:,idxs_beta)); 
title('Beta'); xlabel('Connectivity (Icoh)');

frac_sig_delta = length(find(conspec_sig(:,:,idxs_delta)))/length(find(conspec(:,:,idxs_delta)));
frac_sig_theta = length(find(conspec_sig(:,:,idxs_theta)))/length(find(conspec(:,:,idxs_theta)));
frac_sig_alpha = length(find(conspec_sig(:,:,idxs_alpha)))/length(find(conspec(:,:,idxs_alpha)));
frac_sig_beta = length(find(conspec_sig(:,:,idxs_beta)))/length(find(conspec(:,:,idxs_beta)));

figure('Name', 'Distribution of ICoh Sig. values');
subplot(2, 2, 1); histogram(find(conspec_sig(:,:,idxs_delta)));
title(strcat('Delta (%Sig. =', " ", num2str(frac_sig_delta),')')); 
xlabel('Connectivity (Icoh)');
subplot(2, 2, 2); histogram(find(conspec_sig(:,:,idxs_theta))); 
title('Theta'); title(strcat('Theta (%Sig. =', " ", num2str(frac_sig_theta),')')); 
xlabel('Connectivity (Icoh)');
subplot(2, 2, 3); histogram(find(conspec_sig(:,:,idxs_alpha)));
title('Alpha'); title(strcat('Alpha (%Sig. =', " ", num2str(frac_sig_alpha),')')); 
xlabel('Connectivity (Icoh)');
subplot(2, 2, 4); histogram(find(conspec_sig(:,:,idxs_beta))); 
title('Beta'); title(strcat('Beta (%Sig. =', " ", num2str(frac_sig_beta),')')); 
xlabel('Connectivity (Icoh)');