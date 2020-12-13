function features_avg = average_frequency(features, f_vector, bands)

n_bands = size(bands, 2);
siz = size(features);
siz_avg = [n_bands siz(2:end)];
features_avg = zeros(siz_avg);

for b = 1 : n_bands

    % Assign frequency bands indices 
    idxs_band = logical(f_vector >= bands(1,b) ...
        & f_vector < bands(2,b));

    % Compute average power across current frequency band 
    features_avg(b, :, :, :) = mean(features(idxs_band, ...
        :, :, :), 1);

end

% Put bands as the last dimension 
features_avg = shiftdim(features_avg, 1);

end