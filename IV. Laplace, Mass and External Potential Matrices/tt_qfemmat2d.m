% 
% Generate the QTT of Laplace and Mass Matrices in 2-D FEM
% 
% QTT = TT_QFEMMAT2D(D, VARARGIN)
%   Return the exact QTT representation of Laplace matrix, mass matrix and 
%   sum of them
% 
%   [Input Argument]
%       d - Scalar, level of the mesh grid size, total element number will
%           be (2^d+1)^2
%       dom_size - Vector, 2 x 1, optional, size of the domain, i.e. range
%                  of X and Y (default: [1 1])
%       mat_type - Cell or scalar, optional, should be a subset of 
%                  {'Laplace', 'Mass', 'Laplace+Mass'} that are matrix 
%                  types (default: {'Laplace', 'Mass', 'Laplace+Mass'})
% 
%   [Ouput Argument]
%       A, M, K - TT-Matrix, optional, QTT of the Laplace matrix, mass 
%                 matrix and the sum of them resp., outputs will be in the 
%                 order as mat_type
% 
% Details:
%   1. A, M, K have QTT decomposition as below:
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
%                      [ cJ'+bI+cJ ]             [ fJ'+eI+fJ ]
%               H_0 = a[     cJ    ]      H_1 = a[     fJ    ]
%                      [     cJ'   ]             [     fJ'   ]
% 
%                 I = [1 0]                 J = [0 1]
%                     [0 1]                     [0 0]
% 
%      O is arbitary zero matrix，
%      for A, 
% 
%           a = 1 / (3 * hx * hy)
%           b = 4 * (hx^2 + hy^2)
%           c = hy^2 - 2 * hx^2
%           e = hx^2 - 2 * hy^2
%           f = -(hx^2 + hy^2) / 2
% 
%      for M, 
% 
%           a = hx * hy / 36
%           b = 16
%           c = 4
%           e = 4
%           f = 1
% 
%      for K, 
% 
%           a = 1
%           b = 4 * (hx^2 + hy^2) / (3 * hx * hy) + 4 * hx * hy / 9
%           c = (hy^2 - 2 * hx^2) / (3 * hx * hy) + hx * hy / 9
%           e = (hx^2 - 2 * hy^2) / (3 * hx * hy) + hx * hy / 9
%           f = -(hx^2 + hy^2)/ (6 * hx * hy) + hx * hy / 36
% 
%      where hx/hy is size of each element on X/Y axis


function varargout = tt_qfemmat2d(d, varargin)
    % Input number check
    switch nargin
        case 1
            dom_size = [1 1];
            mat_type = {'Laplace', 'Mass', 'Laplace+Mass'};
        case 2
            if ischar(varargin{1})
                dom_size = [1 1];
                mat_type = varargin{1};
            else
                dom_size = varargin{1};
                mat_type = {'Laplace', 'Mass', 'Laplace+Mass'};
            end
        case 3
            dom_size = varargin{1};
            mat_type = varargin{2};
        otherwise
            error('Number of inputs should be 1~3!');
    end


    % Input type check
    if ~(isscalar(d) && isnumeric(d))
        error('Level of the grid size should be a scalar!');
    elseif rem(d, 1) ~= 0
        error('Level of the grid size should be an integer!');
    elseif d <= 0
        error('Level of the grid size should be positive!');

    elseif ~(isvector(dom_size) && isnumeric(dom_size))
        error('Size of domain should be a numeric vector!');
    elseif length(dom_size) ~= 2
        error('Size of domain should be a length-2 vector!');
    elseif any(dom_size <= 0)
        error('Size of domain should be positive!');

    elseif ~iscell(mat_type)
        if ischar(mat_type)
            mat_type = {mat_type};
        else
            error('Types of matrices should be a cell or character!');
        end
    elseif isempty(mat_type) || length(mat_type) > 3
        error('Types of matrices should be non-empty and has length 1~3!');
    else
        for i = mat_type
            if ~ischar(i{:})
                error('Type of matrix should be a character!');
            elseif ~strcmp(i{:}, 'Laplace') && ~strcmp(i{:}, 'Mass') && ...
                   ~strcmp(i{:}, 'Laplace+Mass')
                error(['Type of matrix should be ''Laplace'',', ...
                       ' ''Mass'' or ''Laplace+Mass''']);
            end
        end

    end
    
    
    % Size of element
    hx = dom_size(1)/(2^d+1);
    hy = dom_size(2)/(2^d+1);
    

    % Generate QTT
    for i = mat_type
        % Generate QTT of Laplace matrix
        if strcmp(i{:}, 'Laplace')
            A = tt_qtridiag2d(d, ...
                              [4*(hx^2+hy^2), hy^2-2*hx^2, ...
                               hx^2-2*hy^2, -(hx^2 + hy^2)/2]/(3*hx*hy));

        % Generate QTT of Mass matrix
        elseif strcmp(i{:}, 'Mass')
            M = tt_qtridiag2d(d, [16 4 4 1]*hx*hy/36);

        % Generate QTT of sum of Laplace and mass matrices
        else
            K = tt_qtridiag2d(d, ...
                              [4*(hx^2+hy^2), hy^2-2*hx^2, hx^2-2*hy^2, ...
                               -(hx^2 + hy^2)/2]/(3*hx*hy)+ ...
                              [16 4 4 1]*hx*hy/36);

        end
    end


    % Return outputs
    varargout = cell(length(mat_type), 1);
    for i = 1:length(mat_type)
        switch mat_type{i}
            case 'Laplace'
                varargout{i} = A;
            case 'Mass'
                varargout{i} = M;
            case 'Laplace+Mass'
                varargout{i} = K;
        end
    end


end