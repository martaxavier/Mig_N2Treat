function [Cxy_tri] = vec2tri(Cxy, dim_X, n_X, upper_or_lower)

% Read dimension of the problem 
siz = size(Cxy);

% Create the permutation vector 
perm = 1 : 1 : ndims(Cxy);
perm(dim_X) = [];
perm = [dim_X perm];

% Bring target vector to the first dimension and 
% concatenate the remaining dimensions into a third
siz = siz(perm);
Cxy = reshape(permute(Cxy, perm), [siz(1) prod(siz(2:end))]);

% Find the diagonal indices
[iY,iX] = meshgrid(1 : n_X, 1 : n_X); 
Cxy_tri = zeros(n_X*n_X, prod(siz(2:end)));

switch upper_or_lower 
    % Upper triangular 
    case 'upper'
        Cxy_tri(logical(iX <= iY), :) = Cxy;
    % Lower triangular    
    case 'lower'
        Cxy_tri(logical(iX >= iY), :) = Cxy;
end

% Reshape into a triangular matrix 
Cxy_tri = reshape(Cxy_tri, [n_X, n_X, siz(2:end)]);

% Permute back to the original dimensions 
dim_Y = dim_X + 1; 
perm2 = 1 : 1 : ndims(Cxy_tri);
perm2(1 : dim_X - 1) = dim_Y : perm2(end) - 1;
perm2(dim_X) = 1; perm2(dim_Y) = 2;
Cxy_tri = permute(Cxy_tri, perm2);             

end