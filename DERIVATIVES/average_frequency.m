function features_avg = average_frequency(features, f_vector, bands)

% Number of frequency bands 
n_bands = size(bands,2);

% Dimension of the features matrix
siz = size(features);

% Pre-allocate output matrix 
features_avg = zeros(siz);

% Go through frequency bands 
for b = 1 : n_bands

    % Find frequency indices of the current frequency 
    % band 
    band_idxs = logical(f_vector >= bands(1,b) & ...
        f_vector < bands(2,b));

    % Compute average power across current frequency band 
    features_avg(b, :, :, :) = mean(features(band_idxs, :, :, :), 1);

end

% Permute the output so as to ensure that the frequency 
% bands are the last dimension 
features_avg = permute(features_avg, [2:length(siz) 1]);
            
end