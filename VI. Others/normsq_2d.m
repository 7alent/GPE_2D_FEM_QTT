% 
% Calculate Square of the Norm of a Function
% 
% OUTPUT = NORMSQ_2D(V, INTERVAL, NORM_TYPE, RULE, VARARGIN)
%   Given a 2-D functions in the form of a vector/matrix (function values 
%   on an equidistant grid), or a function expression (a function handle/
%   sym), calulate the square of the L2/H1 norm or a general positive-
%   definite-matrix norm of it
%   Consider Dirichlet boundary condition and suppose the domain of both 
%   functions is xy_bound = [xa, xb] x [ya, yb] 
%   The square of norm is calculated via Gauss quadrature on each element 
%   if at least one of the input functions (v and w) is a vector, while
%   calculation via quadratic form is also optional
%   A simple version of norm_sq2()
% 
%   [Input Argument]
%       v - Vector/matrix or function handle/sym, maybe TT in the future
%       interval - Vector, 4 x 1, the 1-st and 2-nd entries are the
%                  left/right endpoints of the X axis while last 2 entries
%                  denotes the range of Y axis similarly
%       norm_type - Character or matrix, 'L2' or 'H1' for L2/H1 norm resp.,
%                   set it as a positive definite matrix to calculate
%                   matrix norm
%       rule - Vector or cell, rules of Gauss quadrature of x/y, can be 
%              given as an 2 x 1 integer vector (orders of quadrature 
%              w.r.t. x and y) or a 2 x 1 cell (its values are 2-column 
%              matrices, whose columns are abscissas and weights of x/y 
%              resp.)
%       method - Character, optional, which method to use to calculate the 
%                square of the norm, only works when v is a vector/matrix 
%                and norm_type is 'L2'/'H1' or a general positive definite 
%                matrix, 1/0 for 'quadrature by element'/'quadratic form' 
%                method resp.
% 
%   [Ouput Argument]
%       output - Double, square of the norm of the input function
% 
% Details:
%   1. If v is a matrix, it should be the function values on the grid 
%      (boundary points included, i.e. a matrix whose first&last rows and 
%      columns are 0); if it's a vector, it should be the vectorization of 
%      the matrix (by row) we'd just mentioned but without its boundary
%      points (i.e. without first&last rows and columns); if it's a sym, it
%      will be converted to a function handle
%   2. If norm_type is a matrix, method will be switched to 0


function output = normsq_2d(v, interval, norm_type, rule, varargin)
    % Deal with sym input
    if isa(v, 'sym')
        v = matlabFunction(v, 'Vars', {'x','y'});
    end
    

    % Optional input argument
    if nargin == 5
        method = varargin{1};
    elseif nargin == 4
        method = 1;
    else
        error('Number of arguments should be 4 or 5');
    end


    % Norm type check
    if ~ischar(norm_type)
        method = 0;
    elseif ischar(norm_type)
        if ~strcmp(norm_type, 'L2') && ~strcmp(norm_type, 'H1')
            error('If norm_type is a character, it should be L2/H1');
        end
    end

    
    % Only L2/H1 norm is available for function handle
    if isa(v, 'function_handle')
        method = 1;
        if ~ischar(norm_type)
            error('For function handle input, set norm_type = L2/H1');
        end
    end
    

    % Form conversion before calculation
    if method
        % 'Quadrature by elements' method only works for matrices
        if ~isa(v, 'function_handle') && isvector(v)
            v = reshape(v, [sqrt(numel(v)) sqrt(numel(v))])';
            v = [zeros(1, size(v, 2)+2);
                 zeros(size(v, 1), 1), v, zeros(size(v, 1), 1);
                 zeros(1, size(v, 2)+2)];
        end

    else
        % 'Quadratic form' method only works for vectors
        if ~isvector(v)
            v = v(2:end-1, 2:end-1);
            v = reshape(v', [numel(v) 1]);
        end

        % Generate Laplacian/mass matrices if needed
        if ischar(norm_type) && isvector(v)
            if strcmp(norm_type, 'L2')
                norm_type = fem_mat_2d_boost(sqrt(size(v, 1))+1, ...
                                             interval, 'Mass', rule);
            else
                norm_type = fem_mat_2d_boost(sqrt(size(v, 1))+1, ...
                                             interval, 'Laplace+Mass', ...
                                             rule);
            end
        end

    end


    % Calculate the square of the norm
    if isa(v, 'function_handle')
        if strcmp(norm_type, 'L2')
            output = gq2(@(x, y)(v(x, y)^2), interval, rule);
        elseif strcmp(norm_type, 'H1')
            vx = matlabFunction(diff(sym(v), 'x'), 'Vars', {'x','y'});
            vy = matlabFunction(diff(sym(v), 'y'), 'Vars', {'x','y'});
            output = gq2(@(x, y)(v(x, y)^2+vx(x, y)^2+vy(x, y)^2), ...
                         interval, rule);
        end
    
    elseif method % 'Quadrature by element'
        % Some pre-calculation
        N = size(v, 1)-1; % Number of elements on each row
        M = size(v, 2)-1; % Number of elements on each column
        x_range = interval(2)-interval(1);
        y_range = interval(4)-interval(3);
        hx = x_range/N;
        hy = y_range/M;
        J = hx*hy/4; % Jacobian
        
        % Shape function and its gradients
        f = @(pn, qm, p, q)((1+pn*p)*(1+qm*q)/4);
        if strcmp(norm_type, 'H1')
            fp = @(pn, qm, p, q)((pn*(1+qm*q))/4*2/hx);
            fq = @(pn, qm, p, q)(((1+pn*p)*qm)/4*2/hy);
        end

        % Calculate the integral on each element
        output = 0;
        for m = 1:M % Each element on a row
        for n = 1:N % Each element on a column
            r = @(p, q)(v(m, n)*f(-1, -1, p, q) ...
                        +v(m, n+1)*f(1, -1, p, q) ...
                        +v(m+1, n)*f(-1, 1, p, q) ...
                        +v(m+1, n+1)*f(1, 1, p, q));
            rp = @(p, q)(v(m, n)*fp(-1, -1, p, q) ...
                         +v(m, n+1)*fp(1, -1, p, q) ...
                         +v(m+1, n)*fp(-1, 1, p, q) ...
                         +v(m+1, n+1)*fp(1, 1, p, q));
            rq = @(p, q)(v(m, n)*fq(-1, -1, p, q) ...
                         +v(m, n+1)*fq(1, -1, p, q) ...
                         +v(m+1, n)*fq(-1, 1, p, q) ...
                         +v(m+1, n+1)*fq(1, 1, p, q));
            if strcmp(norm_type, 'L2')
                integrand = @(p, q)(r(p, q)^2*J);
                output = output+gq2(integrand, [-1 1 -1 1], rule);
            else
                integrand = @(p, q)((r(p, q)^2+rp(p, q)^2+rq(p, q)^2)*J);
                output = output+gq2(integrand, [-1 1 -1 1], rule);
            end
        end
        end

    else % 'Quadratic form'
        output = v'*norm_type*v;
    end


end