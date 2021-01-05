function [Cxy] = tril2symmetric(Cxy)

% Read dimension of the problem 
[n_pnts, n_chans, ~, n_bands]= size(Cxy);

% Flip the lower triangular matrix 
% to build a symmetric matrix 
Cxy_tril = Cxy;
Cxy_T = permute(Cxy,[1, 3, 2, 4]);
Cxy = Cxy_T + Cxy;

% Keep the diagonal as it was in
% the original Cxy matrix 
Cxy = reshape(permute(Cxy, [3, 2, 1, 4]), ...
    [n_chans n_chans n_pnts*n_bands]);
Cxy_tril = reshape(permute(Cxy_tril, ...
    [3, 2, 1, 4]), [n_chans n_chans n_pnts*n_bands]);

% Extract the diagonal indices in the 
% reshaped matrices, and replace the 
% diagonal elements of the symmetric
% matrices by the diagonal elements of 
% the lower triangular matrices 
diag_idxs = find(repmat(eye(n_chans, ...
    n_chans),[1 1 n_pnts*n_bands]));
Cxy(diag_idxs) = Cxy_tril(diag_idxs);

% Resshape the connectivity matrix into 
% its original size 
Cxy = permute(reshape(Cxy,[n_chans ...
    n_chans n_pnts n_bands]), [3, 2, 1, 4]);