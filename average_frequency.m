function features_avg = average_frequency(features, f_vector, bands, dim_avg)

% Read size of the problem 
n_bands = size(bands, 2);

% Shift dimensions until frequency 
% is the first dimension 
features = shiftdim(features, ...
    dim_avg - 1);
siz_avg = size(features);
siz_avg(1) = n_bands;
features_avg = zeros(siz_avg);

for b = 1 : n_bands

    % Assign frequency bands indices 
    idxs_band = logical(f_vector >= bands(1, b) ...
        & f_vector < bands(2, b));

    % Compute average power across current frequency band 
    features_avg(b, :, :, :) = mean(features(idxs_band, ...
        :, :, :), 1);

end

features_avg = shiftdim(features_avg, ...
    length(siz_avg) - dim_avg + 1);

end