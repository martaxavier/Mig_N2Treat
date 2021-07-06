function [Cxy, p_values, f_vector] = get_connectivity(X, fs_X, max_freq_res, metric)

p_values = [];

switch metric

    % Imaginary part of coherence
    case 'icoh'

        [Cxy, p_values, f_vector, ~, ~, ~] = ...
            bst_cohn(X, X, fs_X, max_freq_res, 0.5, ...
            'icohere', 1, [], 100);
        Cxy = squeeze(Cxy);

    % Weighted phase lag index 
    % (debiased) 
    case 'wpli'

        [Cxy, f_vector, ~, ~, ~] = ...
            wpli(X, X, fs_X, max_freq_res, 0.5, 1);
        Cxy = squeeze(Cxy);

end

end