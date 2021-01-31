function [Cxy] = tri2symmetric(Cxy_tri, dim_X, dim_Y)

% Read dimension of the problem 
siz = size(Cxy_tri);
n_X = siz(dim_X);

% Create the permutation vector 
perm1 = 1 : 1 : ndims(Cxy_tri);
perm1(dim_X) = dim_Y;
perm1(dim_Y) = dim_X;

% Flip the triangular matrix 
% to build a symmetric matrix 
Cxy_tri_T = permute(Cxy_tri, perm1);
Cxy = Cxy_tri + Cxy_tri_T;

% Put the symmetrical matrix in the 
% first two dimensions and concatenate
% the remaining dimensions into a third 
perm1([dim_X dim_Y]) = []; 
perm2 = [dim_X dim_Y perm1];

% Keep the diagonal as it was in the original Cxy matrix 
siz = siz(perm2);
Cxy = reshape(permute(Cxy, perm2), [n_X n_X prod(siz(3:end))]);
Cxy_tri = reshape(permute(Cxy_tri, perm2), [n_X n_X prod(siz(3:end))]);

% Find the diagonal indices and change them 
[iY, iX, ~] = ndgrid(1 : n_X, 1 : n_X, 1 : prod(siz(3:end))); 
Cxy(logical(iX == iY)) = Cxy_tri(logical(iX == iY));

% Reshape and permute back to the original shape 
[~, perm3] = sort(perm2);
Cxy = permute(reshape(Cxy, [n_X, n_X, siz(3:end)]), perm3);

end