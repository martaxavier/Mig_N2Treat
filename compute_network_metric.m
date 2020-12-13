function [Gxy] = compute_network_metric(Cxy, metric)

% Read dimension of the problem 
[n_bands,n_pnts,n_chans,~]= size(Cxy);

% Flip the lower triangular matrix 
% to build a symmetric matrix 
Cxy_tril = Cxy;
Cxy_T = permute(Cxy,[1,2,4,3]);
Cxy = Cxy_T + Cxy;

% Keep the diagonal as it was in
% the original Cxy matrix 
Cxy = reshape(permute(Cxy,[4,3,1,2]), ...
    [n_chans n_chans n_bands*n_pnts]);
Cxy_tril = reshape(permute(Cxy_tril, ...
    [4,3,1,2]),[n_chans n_chans n_bands*n_pnts]);

% Extract the diagonal indices in the 
% reshaped matrices, and replace the 
% diagonal elements of the symmetric
% matrices by the diagonal elements of 
% the lower triangular matrices 
diag_idxs = find(repmat(eye(n_chans, ...
    n_chans),[1 1 n_bands*n_pnts]));
Cxy(diag_idxs) = Cxy_tril(diag_idxs);

% Resshape the connectivity matrix into 
% its original size 
Cxy = permute(reshape(Cxy,[n_chans ...
    n_chans n_bands n_pnts]),[3,4,1,2]);
Gxy = zeros(n_bands, n_pnts, n_chans);

switch metric 

    % Weighted node degree
    case 'wnd'
        
        Gxy = sum(Cxy,4);
     
    % Clustering coefficient     
    case 'cc'
        
        for b = 1 : n_bands
            
            for p = 1: n_pnts
                    
                    Gxy(b, p, :) = ...
                        clustering_coef_wu(Cxy(:, :, b, p));
            end
        
        end

end 

    end
