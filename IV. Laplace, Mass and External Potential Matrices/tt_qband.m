% 
% Calculate the QTT of Band Matrices
% 
% QTT = TT_QBAND(A, TOL)
%   Get the QTT representation of a given sparse band matrix
% 
%   [Input Argument]
%       A - Matrix, band matrix of size N x N with N = 2^d
%       tol - Double, tolerance of QTT approximation
% 
%   [Ouput Argument]
%       qtt - Cell Array, QTT representation of the band matrix
% 
% Details:
%   1. We first define the 'k-th super-diagonal' (a N x 1 vector):
%
%          a_k = (0, ... , 0, A(1, k+1), A(2, k+2), ... , A(N-k, N))'
%         a_{-k} = (A(k+1, 1), A(k+2, 2), ... , A(N, N-k), 0, ... , 0)'
%
%      here k = 1, ... , N-1 (it coincides with the case k = 0 for the 
%      main diagonal a_0 = (A(1, 1), ... , A(N, N))').
%      Then we define the 'k-th' shift matrix as:
%
%               S_k(i, j) = 1 when j=i+k and S_k(i, j) = 0 o.w.
%
%      here k = -(N-1), ... , -1, 0, 1, ... , N-1.
%      We denote the diagonalization of a vector by diag(.), and let L/U be
%      the lower/upper triangular part of A (without main diagonal a_0) 
%      then we get:
%
%                            A = diag(a_0) + L + U
%         L = S_{-(N-1)} * diag(a_{-(N-1)}) + ... + S_{-1} * diag(a_{-1})
%               U = S_1 * diag(a_1) + ... + S_{N-1} * diag(a_{N-1})
%      
%      specially, if A is symmetric:
%      
%                            A = diag(a_0) + U' + U
%      
%      If we get the QTT approximation of diagonal vectors v_k and exact
%      QTT decomposition of shift matrices S_k, we can use QTT addition, 
%      multiplication and transposition to calculate the QTT approximation 
%      of A
%      This method is reasonable for the case that A is a band matrix and 
%      only a few number of diagonals are non-zero so we won't have to 
%      conduct too many times of QTT operations to get the QTT of A


function qtt = tt_qband(A, tol)
    % Get the level
    N = size(A, 1);
    d = log(N)/log(2);
    

    % Symmetry check
    is_symmetric = issymmetric(A);
    if is_symmetric
        A = triu(A); % Upper triangular part of A
    end


    % Get non-zero diagonals
    [diagonal_vectors, diagonal_indices] = spdiags(A);


    % Get QTT of diagonals and multiply them with shift matrices
    diagonal_cell = cell(length(diagonal_indices), 1);
    for i = 1:length(diagonal_indices)
        k = diagonal_indices(i);
        a_k = tt_qfromfull(diagonal_vectors(:, i), 1, d, tol, 2);
        if k ~= 0 % Sub/super-diagonals
            S_k = tt_matrix(tt_qshift(d, -k));
            diagonal_cell{i} = round(S_k*diag(tt_tensor(a_k)), tol);
        else % Main diagonal
            diagonal_cell{i} = diag(tt_tensor(a_k));
        end
    end


    % Addition of QTT
    if is_symmetric
        for i = 2:length(diagonal_cell)
            if ~exist('qtt', 'var')
                qtt = diagonal_cell{i};
            else
                qtt = qtt+diagonal_cell{i};
            end
        end
        qtt = diagonal_cell{1}+qtt+qtt';
    else
        for i = 1:length(diagonal_cell)
            if ~exist('qtt', 'var')
                qtt = diagonal_cell{i};
            else
                qtt = qtt+diagonal_cell{i};
            end
        end
    end
    qtt = round(qtt, tol);


end