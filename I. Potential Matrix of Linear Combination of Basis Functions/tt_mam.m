% 
% Multi-array Multiplication in TT-tensor Format
% 
% C = TT_MAM(A, B, VARARGIN)
%   Return the multiplication (contraction) product C of two arbitary 
%   arrays A and B in TT-tensor format:
% 
%                               C = A * B
%   
%   by contraction, B is inserted into one of the mode (i.e. contraction 
%   mode) of A
% 
%   [Input Argument]
%       A, B - TT-tensor, TT-matrix or cell array, multiplicators, should 
%              have at least one compatible mode
%       A_i, B_i - Scalar, optional, modes of A and B to be contracted; if 
%                  both not given, they will be set as the first mode of A
%                  that satisfies the contraction requirement and the first 
%                  mode that meets with the former; if only one of them is 
%                  given (as an integer), A_i and B_i will be the same as 
%                  this integer
%       tol - Scalar, optional, if given, TT-rounding will be conducted on
%             the output (default: [])
% 
%   [Ouput Argument]
%       C - Cell array, the product
% 
% Details:
%   1. If A and B have sizes as
% 
%                          A ∈ R^{a_1, ... , a_m}
%                          B ∈ R^{b_1, ... , b_n}
% 
%      The output will be of sizes as (suppose i-th mode of A and j-th 
%      mode of B is contracted)
% 
%                       C ∈ R^{a_1, ... , a_{i-1}, 
%                              b_1, ... , b_{j-1}, 
%                              b_{j+1}, ... , b_n, 
%                              a_{i+1}, ... , a_m}
% 
%      To calculate the output (TT), for d-th core of C we have
% 
%        C_d(x_d_1, ... , x_d_{i-1}, 
%            y_d_1, ... , y_d_{j-1}, 
%            y_d_{j+1}, ... , y_d_n, 
%            y_d_{i+1}, ... , y_d_m) 
%        = Σ_{z_d}{
%                  A_d(x_d_1, ... , x_d_{i-1}, z_d, x_d_{i+1}, ... , x_d_m)
%                  ⊗
%                  B_d(y_d_1, ... , y_d_{j-1}, z_d, y_d_{j+1}, ... , y_d_n)
%                 }
% 
%      with
% 
%           x_1, ... , x_m are indices w.r.t. modes of d-th core of A
%           y_1, ... , y_n are indices w.r.t. modes of d-th core of B
%           z_d ∈ {1, ... , n_d}, n_d = a_i_d = b_j_d
%           ⊗ is matrix Kronecker product


function C = tt_mam(A, B, varargin)
    % Input check
    if nargin ~= 2 && nargin ~= 3 && nargin ~= 4 && nargin ~= 5
        error('Input number should be 2, 3, 4 or 5!');
    end
    

    % Input type check
    if iscell(A)
        for d = 1:length(A)
            if ~is_array(A{d})
                error(['Inputs should be TT-tensors, TT-matrices or ', ...
                       'cell arrays!']);
            end
        end
    elseif isa(A, 'tt_tensor') || isa(A, 'tt_matrix')
        A = core(A);
    else
        error(['Inputs should be TT-tensors, TT-matrices or ', ...
               'cell arrays!']);
    end

    if iscell(B)
        for d = 1:length(B)
            if ~is_array(B{d})
                error(['Inputs should be TT-tensors, TT-matrices or ', ...
                       'cell arrays!']);
            end
        end
    elseif isa(B, 'tt_tensor') || isa(B, 'tt_matrix')
        B = core(B);
    else
        error(['Inputs should be TT-tensors, TT-matrices or ', ...
               'cell arrays!']);
    end


    % Sizes of A and B
    A_sz = tt_sz(A);
    B_sz = tt_sz(B);


    % Dimensions (number of cores) of A and B
    A_d = length(A);
    B_d = length(B);


    % Dimension check
    if A_d ~= B_d
        error('Incompatible dimensions!');
    end


    % Number of modes of A and B
    A_m = size(A_sz, 2)-2;
    B_m = size(B_sz, 2)-2;


    % Input assignment
    if nargin == 2
        [A_i, B_i] = deal([]);
        tol = [];
    elseif nargin == 3
        if is_array(varargin{1}) && isscalar(varargin{1})
            if mod(varargin{1}, 1) == 0
                [A_i, B_i] = deal(varargin{1});
                tol = [];
            else
                [A_i, B_i] = deal([]);
                tol = varargin{1};
            end
        else
            error('Wrong optional input!');
        end
    else
        A_i = varargin{1};
        B_i = varargin{2};
        tol = varargin{3};
    end


    % Check contraction mode index if it's given
    if ~isempty(A_i) && ~isempty(B_i)
        if ~isscalar(A_i) || ~isnumeric(A_i) || ... 
           ~isscalar(B_i) || ~isnumeric(B_i)
                error('Contraction indices should be numeric scalars!');
        elseif mod(A_i, 1) ~= 0 || A_i < 1 || mod(B_i, 1) ~= 0 || B_i < 1
                error('Contraction indices should be positive integers!');
        end
    end


    % Determine contraction mode index
    n = zeros(A_d, 1); % Flag variable, sizes of contration modes
    if isempty(A_i) && isempty(B_i)
        for A_i = 1:A_m
            for B_i = 1:B_m
                if all(A_sz(:, A_i) == B_sz(:, B_i))
                    n = A_sz(:, A_i);
                    break
                end
            end
            if all(n ~= 0)
                break
            end
        end
    elseif A_i > A_m || B_i > B_m
        error('Contraction indices should be mode indices!');
    elseif A_sz(:, A_i) == B_sz(:, B_i)
        n = A_sz(A_i);
    end


    % Contraction mode size check
    if any(n == 0)
        error('Incompatible contraction modes!');
    end


    % Check rounding tolerance if it's given
    if ~isempty(tol)
        if ~is_array(tol) || ~isscalar(tol)
            error('Tolerance should be a numeric scalar!');
        elseif tol <= 0
            error('Tolerance should be positive!');
        end
    end

    
    % Initialization of product
    C = cell(A_d, 1);
    for d = 1:A_d
        C{d} = zeros([A_sz(d, 1:(A_i-1)) B_sz(d, 1:B_m ~= B_i) ...
                      A_sz(d, (A_i+1):A_m) ...
                      A_sz(d, A_m+(1:2)).*B_sz(d, B_m+(1:2))]);
    end


    % Calculation of product
    for d = 1:A_d
        % Permute and reshape cores
        a = A{d};
        b = B{d};
        a_r = A_sz(d, A_m+(1:2));
        b_r = B_sz(d, B_m+(1:2));
        a = permute(a, [1:(A_i-1) (A_i+1):A_m A_i A_m+(1:2)]);
        a = reshape(a, [numel(a)/(A_sz(d, A_i)*prod(a_r)) ...
                        A_sz(d, A_i) a_r]);
        b = permute(b, [B_i 1:(B_i-1) (B_i+1):B_m B_m+(1:2)]);
        b = reshape(b, [B_sz(d, B_i) numel(b)/(B_sz(d, B_i)*prod(b_r)) ...
                        b_r]);

        % Kronecker product by block
        tmp = zeros([size(a, 1) size(b, 2) a_r.*b_r]);
        for i = 1:size(a, 1)
        for j = 1:size(b, 2)
        for k = 1:size(a, 2)
            tmp(i, j, :, :) = ...
                tmp(i, j, :, :)+...
                reshape(kron(reshape(a(i, k, :, :), a_r), ...
                             reshape(b(k, j, :, :), b_r)), ...
                        [1 1 a_r.*b_r]);
        end
        end
        end
        tmp = reshape(tmp, [A_sz(d, 1:A_m ~= A_i) B_sz(d, 1:B_m ~= B_i) ...
                            a_r.*b_r]);
        tmp = permute(tmp, [1:(A_i-1) A_m:(A_m+B_m-2) A_i:(A_m-1) ...
                            A_m+B_m-1 A_m+B_m]);

        % Assemble the product
        C{d} = tmp;
    end
    C{1} = squeeze(C{1});


    % Rounding
    if ~isempty(tol)
        C = tt_round(C, tol);
    end


end

