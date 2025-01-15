% 
% Generate Multi-diagonal Matrices
% 
% MAT = MULTIDIAG(X, N)
%   Generate the n x n multi-diagonal matrix by given a generator vector x
% 
%   [Input Argument]
%       x - Vector, should have odd length m <= n
%       n - Integer, size of the matrix
% 
%   [Ouput Argument]
%       mat - Matrix, the output multi-diagonal matrix
% 
% Details:
%   1. Given a x = (x_1, ... , x_{k-1}, x_k, x_{k+1}, ... , x_m)' where 
%   k-1 = (m-1)/2, it means the subdiagonals indexed by -(k-1), ... , -1 
%   are filled with x_1, ... , x_(k-1), similarly, superdiagonals indexed 
%   by 1, ..., k-1 are filled with x_{k+1}, ... , x_m, while the diagonal 
%   indexed by 0 is filled with x_k


function mat = multidiag(x, n)
    mat = sparse(n, n);
    k = (length(x)-1)/2+1;
    for i = 1:length(x)
        if i == k % Diagonal
            mat = mat+spdiags(repmat(x(i), [1 n]), 0, n, n);
        else
            mat = mat+spdiags(repmat(x(i), [1 n-1]), i-k, n, n);
        end
    end
end