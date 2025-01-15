% 
% Generate the QTT of Tridiagonal Block Matrices in 2-D FEM
% 
% TT = TT_QTRIDIAG2D(D, X)
%   Return the exact QTT representation of a tridiagonal block matrix whose
%   blocks are also tridiagonal matrices, and all main blocks are the same,
%   all subdiagonal and superdiagonal blocks are the same as well
%   That is
% 
%                       A = tridiag(A_1, A_0, A_1)
%                      A_0 = tridiag(x_2, x_1, x_2)
%                      A_1 = tridiag(x_4, x_3, x_4)
% 
%   [Input Argument]
%       d - Scalar, level of the mesh grid size, total element number will
%           be (2^d+1)^2
%       x - Scalar, x = [x_1 x_2 x_3 x_4] specifes main and sub/super
%           diagonal entries of main and sub/super diagonal blocks
% 
%   [Ouput Argument]
%       TT - TT-Matrix, QTT of the output matrix
% 
% Details:
%   1. A has QTT decomposition as below:
% 
%           F * G^(d-2) * [H_0 H_1] * [I O] * [I O]^(d-2) * [I]
%                                     [O F]   [O G]         [Q]
% 
%      here * denotes rank core product (⋈) and ^ denotes its power, and 
% 
%                 F = [ I  J' J ]
% 
%                     [ I  J' J ]               [ J'+J  ]
%                 G = [ O  J  O ]           Q = [   J   ]
%                     [ O  O  J']               [   J'  ]
% 
%                      [ x_2*(J'+J)+x_1*I ]
%               H_0 =  [       x_2*J      ]
%                      [       x_2J'      ]
% 
%                      [ x_4*(J'+J)+x_3*I ]
%               H_1 =  [       x_4*J      ]
%                      [       x_4J'      ]
% 
%                 I = [1 0]                 J = [0 1]
%                     [0 1]                     [0 0]
% 
%      O is an arbitary zero matrix


function tt = tt_qtridiag2d(d, x)
    tt = cell(2*d, 1);
    I = eye(2);
    J = [0 1; 0 0];

    for i = 1:2*d
        if i == 1
            tt{i} =zeros(2, 2, 4);
            tt{i}(:, :, 1) = I;
            tt{i}(:, :, 2) = J'+J;
            tt{i}(:, :, 3) = J;
            tt{i}(:, :, 4) = J';
        elseif 2 <= i && i <= d-1
            tt{i} = zeros(2, 2, 4, 4);
            tt{i}(:, :, 1, 1) = I;
            tt{i}(:, :, 2, 2) = I;
            tt{i}(:, :, 3, 2) = J';
            tt{i}(:, :, 4, 2) = J;
            tt{i}(:, :, 3, 3) = J;
            tt{i}(:, :, 4, 4) = J';
        elseif i == d
            tt{i} = zeros(2, 2, 4, 2);
            tt{i}(:, :, 1, 1) = I;
            tt{i}(:, :, 2, 2) = I;
            tt{i}(:, :, 3, 2) = J';
            tt{i}(:, :, 4, 2) = J;
        elseif i == d+1
            tt{i} = zeros(2, 2, 2, 3);
            tt{i}(:, :, 1, 1) = x(3)*(J'+J)+x(1)*I;
            tt{i}(:, :, 1, 2) = x(3)*J;
            tt{i}(:, :, 1, 3) = x(3)*J';
            tt{i}(:, :, 2, 1) = x(4)*(J'+J)+x(2)*I;
            tt{i}(:, :, 2, 2) = x(4)*J;
            tt{i}(:, :, 2, 3) = x(4)*J';
        elseif d+2 <= i && i <= 2*d-1
            tt{i} = zeros(2, 2, 3, 3);
            tt{i}(:, :, 1, 1) = I;
            tt{i}(:, :, 2, 1) = J';
            tt{i}(:, :, 3, 1) = J;
            tt{i}(:, :, 2, 2) = J;
            tt{i}(:, :, 3, 3) = J';
        else % i == 2*d
            tt{i} = zeros(2, 2, 3);
            tt{i}(:, :, 1) = I;
            tt{i}(:, :, 2) = J';
            tt{i}(:, :, 3) = J;
        end
    end

    tt = tt_matrix(tt);
end