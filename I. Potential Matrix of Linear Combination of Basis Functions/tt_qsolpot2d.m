% 
% QTT of Stiffness Matrix of Linear Combination of bases in 2-D FEM
% 
% TTM = TT_QSOLPOT2D(V, XY_EP, VARARGIN)
%   Compute the QTT representation of the stiffness matrix of a potential 
%   function v which is a linear combination of basis fucntions in 1-D FEM 
%   (e.g. the FEM solution)
%   That is
% 
%                          v = Σ(v_i * f_i(x, y))
% 
%      f_i(x, y) = [1-|x-x_i|/h_x][1-|y-y_i|/h_y] (i∈{1, ... , 2^(2*d)}) is
%      the basis function w.r.t. node i (in global indexing)
% 
%      h_x, h_y are spacing of each element on X, Y axis
% 
%      Ω = [xa, xb] x [ya, yb] is the domain
% 
%   [Input Argument]
%       v - TT-tensor, cell arrary, vector, function handle or sym, if it's
%           a (TT-)vector, it should be length-2^d; if v is a function 
%           handle or sym, it will be approximated by linear combination of
%           bases of the grid and the coefficient vector (e.g. linear 
%           combination coefficient at each node) will be approximated by a
%           QTT via TT-SVD
%       xy_ep - Vector, optional, set it as [xa xb ya yb] if the domain is 
%               [xa, xb] x [ya, yb]
%       d - Scalar, optional, not less than 2, level of the grid, vector v
%           will be length-2^(2*d)
%       tol - Scalar, optional, tolerance of calculating the QTT of vector 
%             v (default: 1e-14)
% 
%   [Ouput Argument]
%       ttm - TT-matrix, QTT of the output stiffness matrix B
% 
% Details:
%   1. The stiffness matrix B is computed by
% 
%                              B = A * v
% 
%       where
% 
%                    A = sepdiag(a, b, b, b, b, b, b)
% 
%       with
% 
%           sepdiag() denotes septuple diagonal 3 level tensor
% 
%           a = h_x / 2
% 
%           b = 5 * h_x / 12
% 
%       Since A is 'perfect shuffle' symmetric, it's first reshaped into a
%       matrix then multiplied by v


function ttm = tt_qsolpot2d(v, xy_ep, varargin)
    % Input number check
    if nargin ~= 2 && nargin ~= 3 && nargin ~= 4
        error('Number of inputs should be 2, 3 or 4!');
    end
    

    % Input initialization
    [d, tol, h_x, h_y] = deal([]);


    % Input assignment
    if nargin > 2
        for i = 1:length(varargin)
            if is_array(varargin{i}) && isscalar(varargin{i})
                if varargin{i} <= 0
                    error('Optional inputs should be positive!');
                elseif mod(varargin{i}, 1) == 0
                    if isempty(d)
                        d = varargin{i};
                    else
                        tol = varargin{i};
                    end
                elseif isempty(tol)
                    tol = varargin{i};
                else
                    error('Level of grid should be a positive integer!');
                end
            else
                error('Wrong type of optional inputs!');
            end
        end
    end


    % Input check of endpoints
    if ~isnumeric(xy_ep) || ~isvector(xy_ep)
        error('Endpoints of the interval should be a numeric vector!');
    elseif length(xy_ep) ~= 4
        error('Endpoints of the interval should be length-4!');
    elseif xy_ep(1) >= xy_ep(2) || xy_ep(3) >= xy_ep(4)
            error('Left endpoint should be less than the right one!');
    end
    
    
    % Input type check of level of grid
    if ~isempty(d)
        if ~isnumeric(d) || ~isscalar(d)
            error('Level of grid should be a numeric scalar!');
        elseif mod(d, 1) ~= 0 || d < 2
            error('Level of grid should be an integer not less than 2!');
        elseif ~isnumeric(tol) || ~isscalar(tol)
            error('Tolerance should be a numeric scalar!');
        elseif d <= 0
            error('Tolerance be larger than 0!');
        end
    end
    

    % Input type check, conversion and generation of v
    if isa(v, 'sym')
        v = matlabFunction(v);
    end

    if isa(v, 'function_handle')
        if isempty(d)
            error(['When inputting a function, ', ...
                   'level of grid should be given!']);
        elseif d < 2
            error('Level of grid should be not less than 2!');
        end
        h_x = (xy_ep(2)-xy_ep(1))/(2^d+1);
        h_y = (xy_ep(4)-xy_ep(3))/(2^d+1);
        [x, y] = meshgrid(h_x:h_x:(2^d*h_x), h_y:h_y:(2^d*h_y));
        v = v(x, y);
        v = reshape(v', [numel(v) 1]);
    end

    if iscell(v)
        v = tt_tensor(v);
        n = prod(v.n);
        d = log(n)/log(2)/2;
        if d < 2
            error('Level of grid should be not less than 2!');
        end
        if any(v.n ~= 2)
            error('Only QTT is accpeted!');
        end
    elseif isvector(v)
        d = log(length(v))/log(2)/2;
        if mod(d, 1) ~= 0
            error('Length of vector should be 2^d!');
        elseif d < 2
            error('Level of grid should be not less than 2!');
        end
        v = tt_qfromfull(v, 1, 2*d, tol, 2);
        v = tt_tensor(v);
    else
        d = v.d/2;
        if d < 2
            error('Level of grid should be not less than 2!');
        end
        if any(v.n ~= 2)
            error('Only QTT is accpeted!');
        end
    end
    

    % Spacings of grid
    if isempty(h_x)
        h_x = (xy_ep(2)-xy_ep(1))/(2^d+1);
    end
    if isempty(h_y)
        h_y = (xy_ep(4)-xy_ep(3))/(2^d+1);
    end


    % Generate Tri-basis-integral Tensor
    A = tt_qsepdiag2d(d, [36 2 6 1]*h_x*h_y/144);


    % Compute stiffness matrix
    ttm = tt_mam(A, v, tol);
    

end