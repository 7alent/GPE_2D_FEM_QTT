% 
% Kronecker Product of 3D Tensors
% 
% C = TKRON(A, B)
%   Calculate the Kronecker product C of twon 3 level tensors A and B
%   If A ∈ R^{a_1 x a_2 x a_3} and B ∈ R^{b_1 x b_2 x b_3}, the output will
%   be 
% 
%      C = (A_{ijk} * B) ∈ R^{(a_1 * b_1) x (a_2 * b_2) x (a_3 * b_3)}
% 
%   It's a natrual extensiton of matrix Kronecker product
% 
%   [Input Argument]
%       A, B - Array, 3 level tensors to be operated
% 
%   [Ouput Argument]
%       C - Array, 3 level tensor, Kronecker product of inputs
% 
% Details:
%   1. This function only works for 3 level tensor (array) rather than
%      vectors and matrices


function C = tkron(A, B)
    % Input check
    if ~is_array(A) || ~is_array(B)
        error('Inputs should be arrays!');
    elseif ndims(A) ~= 3 || ndims(B) ~= 3
        error('Inputs should be level-3!');
    end


    % Intialization
    A_sz = size(A);
    B_sz = size(B);
    C = zeros(A_sz.*B_sz);


    % Assemble the output
    for i = 1:A_sz(1)
        for j = 1:A_sz(2)
            for k = 1:A_sz(3)
                C(((i-1)*B_sz(1)+1):(i*B_sz(1)), ...
                  ((j-1)*B_sz(2)+1):(j*B_sz(2)), ...
                  ((k-1)*B_sz(3)+1):(k*B_sz(3))) = A(i, j, k)*B;
            end
        end
    end


end