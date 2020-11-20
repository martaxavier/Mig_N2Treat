function [conspec, pvalue] = ...
    compute_connectivity_metric(cspec, metric, win_samples)

%   Computes the connectivity metric specified in 'metric', from the 
%   input cross-spectrum, 'cross_spectrum'
%
%   INPUTS:
%
%           cross_spectrum  - Cross-spectrum lower triangular matrix 
%           metric          - Connectivity metric 
% 
%   OUTPUTS:
%
%           conspec         - Output connectivity spectrum 
%           pvalue          - p-value for the significance of the 
%                             connectivity spectrum 
% 

% Define the func connectivity and graph metrics 
Cxy_metric = extractBefore(metric,'_');
Gxy_metric = extractAfter(metric,'_');
if ismissing(Cxy_metric); Cxy_metric = metric; end

% Compute the functional connectivity 
[Cxy, pvalue] = get_Cxy(cspec, Cxy_metric, win_samples);

if ismissing(Gxy_metric)
    conspec = Cxy;
    return
end

% Compute the graph metric
conspec = get_Gxy(Cxy, Gxy_metric);

% if this is ok, then use another function for the pvalue? lets see 

end


    function [Cxy, pval] = get_Cxy(cspec, metric, win_samples)

    % Read dimension of the problem 
    [~,~,n_chans,~]= size(cspec);

    % Pre-allocate Cxy
    Cxy = zeros(size(cspec));
    pval = Cxy;

    % Go through channels 
    for c1 = 1 : n_chans

        % Lower triangular matrix 
        for c2 = 1 : c1

            switch metric 

                 % == Imaginary Coherence 
                case 'icoh'

                    % Cross-spectrum of the current pair of channels
                    Sxy = cspec(:, :, c1, c2);
                    Sxx = cspec(:, :, c1, c1);
                    Syy = cspec(:, :, c2, c2);

                    % Coherence - Cxy = Sxy/sqrt(Sxx*Syy)
                    Cxy_chan = bsxfun(@rdivide, Sxy, sqrt(Sxx));
                    Cxy_chan = bsxfun(@rdivide, Cxy_chan, sqrt(Syy));

                    % Imaginary Part of Coherency - ICxy = abs(im(Cxy))
                    Cxy_chan =  abs(imag(Cxy_chan));

                    % P-value for the current pair of channels
                    n_fft = max(256,2^log2(win_samples));
                    pval_chan = max(0, 1 - abs(Cxy_chan).^2) ...
                        .^ floor(win_samples / n_fft);

                % == Phase Locking Value 
                case 'plv'


            end

            Cxy(:, :, c1, c2) = Cxy_chan;
            pval(:, :, c1, c2) = pval_chan;

        end 

    end

    end

    function [Gxy] = get_Gxy(Cxy,metric)

    % Return if there is no 
    % graph metric to compute 
    if isempty(metric)
        return
    end

    % Read dimension of the problem 
    [n_freq,n_pnts,n_chans,~]= size(Cxy);

    % Flip the lower triangular matrix 
    % to build a symmetric matrix 
    Cxy_tril = Cxy;
    Cxy_T = permute(Cxy,[1,2,4,3]);
    Cxy = Cxy_T + Cxy;

    % Keep the diagonal as it was in
    % the original Cxy matrix 
    Cxy = reshape(permute(Cxy,[4,3,1,2]), ...
        [n_chans n_chans n_freq*n_pnts]);
    Cxy_tril = reshape(permute(Cxy_tril, ...
        [4,3,1,2]),[n_chans n_chans n_freq*n_pnts]);

    % Extract the diagonal indices in the 
    % reshaped matrices, and replace the 
    % diagonal elements of the symmetric
    % matrices by the diagonal elements of 
    % the lower triangular matrices 
    diag_idxs = find(repmat(eye(n_chans, ...
        n_chans),[1 1 n_freq*n_pnts]));
    Cxy(diag_idxs) = Cxy_tril(diag_idxs);

    % Resshape the connectivity matrix into 
    % its original size 
    Cxy = permute(reshape(Cxy,[n_chans ...
        n_chans n_freq n_pnts]),[3,4,1,2]);

    switch metric 

        case 'wnd'

            Gxy = mean(Cxy,4);

    end 

    end
