% 
% 2-D Finite Element Laplace and Mass Matrices
% 
% VARARGOUT = FEM_MAT_2D_EXACT(N, VARARGIN)
%   Return Laplace and mass stiffness matrices in 2-D finite element 
%   method (according to their exact expression)
% 
%   [Input Argument]
%       N - Integer, number of elements on X/Y axis, total number of
%           elements will be N^2
%       dom_size - Vector, 2 x 1, optional, size of the domain, i.e. range
%                  of X and Y (default: [1 1])
%       mat_type - Cell or scalar, optional, should be a subset of 
%                  {'Laplace', 'Mass', 'Laplace+Mass'} that are matrix 
%                  types (default: {'Laplace', 'Mass', 'Laplace+Mass'})
% 
%   [Ouput Argument]
%       A, M, K - Matrix, optional, the Laplace matrix, mass matrix and the
%                 sum of them resp., outputs will be in the order as 
%                 mat_type
% 
% Details:
%   1. This code runs faster than fem_mat_2d_boost() because it doesn't 
%      need to assemble the global stiffness matrices using element 
%      sitffness matrices; but for potential stiffness matrix, please use 
%      fem_mat_2d_boost() instead
%   2. If length of each element on X/Y axis is h_x/h_y, then both Laplace
%      matrix A and mass matrix M are blockwise-tridiagonal matrices like:
% 
%                         A = tridiag(A_1, A_0, A_1)
%                         M = tridiag(M_1, M_0, M_1)
% 
%      here A_0, A_1, M_0, M_1 are (N-1) x (N-1) matrices (i.e. blocks) and
%      they are also tridiagonal:
% 
%                       A_0 = a * tridiag(c_0, b_0, c_0)
%                       A_1 = a * tridiag(c_1, b_1, c_1)
%                       M_0 = d * tridiag(4, 16, 4)
%                       M_1 = d * tridiag(1, 4, 1)
% 
%      with
% 
%                      a = 1 / (3 * h_x * h_y)
%                    b_0 = 4 * [(h_x)^2 + (h_y)^2]
%                    c_0 = (h_x)^2 - 2 * (h_y)^2
%                    b_1 = (h_y)^2 - 2 * (h_x)^2
%                    c_1 = -(1/2) * [(h_x)^2 + (h_y)^2]
%                      d = (1/36) * h_x * h_y


function varargout = fem_mat_2d_exact(N, varargin)
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
    if ~(isscalar(N) && isnumeric(N))
        error('Number of elements on X/Y axis should be a scalar!');
    elseif rem(N, 1) ~= 0
        error('Number of elements on X/Y axis should be an integer!');
    elseif N <= 0
        error('Number of elements on X/Y axis should be positive!');

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
    hx = dom_size(1)/N;
    hy = dom_size(2)/N;


    % Generate matrices
    D = multidiag([1 0 1], N-1);
    I = speye(N-1);
    for i = mat_type
        % Generate Laplace matrix
        if strcmp(i{:}, 'Laplace')
            B = multidiag([(hx^2-2*hy^2)/(3*hx*hy), ...
                           4*(hx^2+hy^2)/(3*hx*hy), ...
                           (hx^2-2*hy^2)/(3*hx*hy)], N-1);
            C = multidiag([-(hx^2+hy^2)/(6*hx*hy), ...
                           (hy^2-2*hx^2)/(3*hx*hy), ...
                           -(hx^2+hy^2)/(6*hx*hy)], N-1);
            A = kron(I, B)+kron(D, C);

        % Generate Mass matrix
        elseif strcmp(i{:}, 'Mass')
            B = multidiag([hx*hy/9 4*hx*hy/9 hx*hy/9], N-1);
            C = multidiag([hx*hy/36 hx*hy/9 hx*hy/36], N-1);
            M = kron(I, B)+kron(D, C);

        % Generate sum of Laplace and mass matrices
        else
            B = multidiag([(hx^2-2*hy^2)/(3*hx*hy)+hx*hy/9, ...
                           4*(hx^2+hy^2)/(3*hx*hy)+4*hx*hy/9, ...
                           (hx^2-2*hy^2)/(3*hx*hy)+hx*hy/9], N-1);
            C = multidiag([-(hx^2+hy^2)/(6*hx*hy)+hx*hy/36, ...
                           (hy^2-2*hx^2)/(3*hx*hy)+hx*hy/9, ...
                           -(hx^2+hy^2)/(6*hx*hy)+hx*hy/36], N-1);
            K = kron(I, B)+kron(D, C);

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