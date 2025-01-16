% 
% Calculate L2 or H1 Error between Two Functions
% 
% ERR = FEM_ERR_2D_TT(V, W, ERR_TYPE, VARARGIN)
%   Given 2 2-D functions in the form of a vector (function values on an 
%   equidistant grid) or a function expression, calulate the L2/H1 error 
%   between them
%   Consider Dirichlet boundary condition and suppose the domain of both 
%   functions is xy_bound = [xa, xb] x [ya, yb] 
%   The error is calculated by Gauss quadrature on each element if at least
%   one of the input functions (v and w) is a vector
%   This function is a updated version of fem_err_2d_simple() that allows
%   TT-tensor inputs and calculation of norm (not just error)
% 
%   [Input Argument]
%       v, w - Vector, function handle/sym or TT-tensor, if both are
%              TT-tensors, they should have same number of cores and mode
%              sizes of corresponding core should also be the same,
%              otherwise they will be converted to vectors
%       err_type - Character or (TT-)matrix, 'L2' or 'H1' for L2/H1 error 
%                  resp., or a positive definite matrix for positive-
%                  definite-matrix norm
%       xy_ep - Vector, optional, 4 x 1 vectors, left and right endpoints 
%               of x and y intervals, set it as [xa xb ya yb] for domain
%               [xa, xb] x [ya, yb] (default: [-1 1 -1 1])
%       gqr - Vector or cell, optional, rules of Gauss quadrature of x/y, 
%             can be given as an 2 x 1 integer vector (orders of quadrature 
%             w.r.t. x and y) or a 2 x 1 cell (its values are 2-column 
%             matrices, whose columns are abscissas and weights of x/y 
%             resp.) (default: [5 5])
%       method - Character, optional, method to calculate the square of the
%                norm, only works when both v and w are vectors, 1/0 for 
%                'quadrature by element'/'quadratic form' methods resp., if
%                err_type is a (TT-)matrix this argument will be reset as 0
%                (default: [])
%       tol - Scalar, optional, tolerance of TT rounding (default: 1e-14)
% 
%   [Ouput Argument]
%       err - Double, L2/H1 error of the 2 input functions
% 
% Details:
%   1. To calculte the L2/H1 norm of v/w, set w/v as 0, but using codes
%      below is recommended
%               norm_v = sqrt(normsq_2d(v, xy_ep, 'L2', gqr));
%   2. Gauss quadrature order (gqr) is recommended to be large (at least 20 
%      for a function handle input, especially when it oscillates severely)
%   3. Sym inputs will be converted to function handles before calculation
%   4. If mass matrix M and its QTT already exist, L2 error can be 
%      calculated by codes below
% 
%           err = fem_err_2d_tt(v, w, M);
%           err_qtt = fem_err_2d_tt(v_qtt, w_qtt, M_qtt, tol);
% 
%     which means argument method isn't always needed as long as other 
%     optional inputs (e.g. tol) is correct (i.e. tol ∈ (0, 1)), because 
%     this function first assigns correct inputs to corresponding arguments
%     then assigns the rest to empty arguments in order
%   5. This function treat v and w as two eigen vectors, so it actually
%      calculates norm of v-w and v+w and select their minimum as the
%      output (that's because if v/w is an eigen vector, so is -v/-w); if
%      only norm of v-w/v+w is needed, set one of the vector input as
%      zeros:
% 
%           err = fem_err_2d_tt(v-w, 0, 'L2'); % L2 norm of v-w
%           err = fem_err_2d_tt(v+w, 0, 'H1'); % H1 norm of v+w
% 
% Future Routines:
%   1. Deal with the case that x and y axes have different step length 
%      (i.e. number of elements on each row ~= number of elements on each 
%      column, N ~= M)


function err = fem_err_2d_tt(v, w, err_type, varargin)
    % Input number check
    if nargin < 3 || nargin > 7
        error('Input number should be 3~7!');
    end

    
    % Intialization of optional inputs
    [xy_ep, gqr, method, tol] = deal([]);

    
    % Input assignment
    if nargin > 3
        for i = 1:length(varargin)
            if iscell(varargin{i})
                gqr = varargin{i};
            elseif is_array(varargin{i}) && isvector(varargin{i})
                if isscalar(varargin{i})
                    if varargin{i} == 0 || varargin{i} == 1
                        method = varargin{i};
                    elseif 0 < varargin{i} && varargin{i} < 1
                        tol = varargin{i};
                    elseif isempty(method)
                        method = varargin{i};
                    elseif isempty(tol)
                        tol = varargin{i};
                    end
                elseif length(varargin{i}) == 2
                    gqr = varargin{i};
                elseif length(varargin{i}) == 4
                    xy_ep = varargin{i};
                else
                    error('Wrong type of input!');
                end
            end
        end
    end


    % Defual inputs
    if isempty(xy_ep)
        xy_ep = [-1 1 -1 1];
    end
    if isempty(gqr)
        gqr = [5 5];
    end
    if isempty(tol)
        tol = 1e-14;
    end


    % Input type check for v and w
    if is_array(v)
        if ~isvector(v)
            error('Multi-dimensional input is unacceptable!');
        elseif isrow(v)
            v = v';
        end
    elseif isa(v, 'sym')
        v = matlabFunction(v, 'Vars', {'x','y'});
    elseif ~isa(v, 'function_handle') && ~isa(v, 'tt_tensor')
        error('Wrong type of input for 1st argument!');
    end

    if is_array(w)
        if ~isvector(w)
            error('Multi-dimensional input is unacceptable!');
        elseif isrow(w)
            w = w';
        end
    elseif isa(w, 'sym')
        w = matlabFunction(w, 'Vars', {'x','y'});
    elseif ~isa(w, 'function_handle') && ~isa(w, 'tt_tensor')
        error('Wrong type of input for 2nd argument!');
    end


    % Input type check for error type
    if ischar(err_type)
        if ~strcmp(err_type, 'L2') && ~strcmp(err_type, 'H1')
            error('Error type should be ''L2'' or ''H1''!');
        end
    elseif ~ismatrix(err_type) && ~isa(err_type, 'tt_matrix')
        error('Error type should be ''L2'' or ''H1''!');
    end


    % Input type check for endpoints
    if ~is_array(xy_ep) || ~isvector(xy_ep)
        error('Wrong type of input for endpoints!');
    elseif length(xy_ep) ~= 4
        error('Endpoints should be a length-4 vector!');
    elseif xy_ep(1) >= xy_ep(2) || xy_ep(3) >= xy_ep(4)
        error('Left endpoint should be less than the right one!');
    end


    % Input type check for quadrture rule
    if iscell(gqr)
        if length(gqr) ~= 2
           error(['Quadrature rule should be a 2 x 1 cell array ', ...
                  'or vector!']);
        elseif ~ismatrix(gqr{1}) || ~ismatrix(gqr{2})
            error(['If quadrature rule is a cell array, its elements ', ...
                   'should be matrcies!']);
        elseif size(gqr{1}, 2) ~= 2 || size(gqr{2}, 2) ~= 2
            error(['If quadrature rule is a cell array, its elements ', ...
                   'should be matrcies with 2 columns!']);
        elseif size(gqr{1}, 1) ~= size(gqr{2}, 1)
            error(['If quadrature rule is a cell array, its elements ', ...
                   'should be matrcies with same number of rows!']);
        end
    elseif is_array(gqr) && isvector(gqr)
        if length(gqr) ~= 2
            error('If quadrature is a vector, it should be length-2!');
        elseif any(mod(gqr, 1) ~= 0) || any(gqr < 1)
            error(['If quadrature is a vector, it should be a ', ...
                   'positive integer vector!']);
        end
    end


    % Input type check for error calculation method
    if ~isempty(method)
        if ~is_array(method) || ~isscalar(method)
            error('Method should be a numeric scalar!');
        elseif method ~= 0 && method ~= 1
            error('Method should be 0 or 1!');
        end
    end
    
    
    % Input type check for TT-rounding tolerance
    if ~is_array(tol) || ~isscalar(tol)
        error('Tolerance should be a numeric scalar!');
    elseif tol <= 0
        error('Tolerance should be positive!');
    end
    

    % Format conversion of v and w
    % Let v be the function handle
    % Also ensure any vector to be a column vector
    % Only conduct TT calculation when v and w are TT with same sizes and
    % err_type is a TT-matrix with compatible sizes with v and w
    if isa(v, 'function_handle') && isa(w, 'function_handle')
        input_type = 'Both are functions';
    elseif isa(v, 'function_handle') && ~isa(w, 'function_handle')
        input_type = 'Function and vector';
        w = full(w);
    elseif ~isa(v, 'function_handle') && isa(w, 'function_handle')
        input_type = 'Function and vector';
        [v, w] = deal(w, full(v));
    else
        input_type = 'Both are vectors';
        if isvector(v) || isvector(w)
            v = full(v);
            w = full(w);
        elseif prod(tt_sz(v, 'mode')) ~= prod(tt_sz(w, 'mode')) || ...
               size(tt_sz(v), 2) ~= size(tt_sz(w), 2) || v.d ~= w.d || ...
               ~isa(err_type, 'tt_matrix')
            v = full(v);
            w = full(w);
        elseif any(any(tt_sz(v, 'mode') ~= tt_sz(w, 'mode'))) || ...
               any(err_type.n ~= err_type.m) || v.d ~= err_type.d
            error('Incompatible sizes of vectors and matrix!');
        elseif any(v.n ~= err_type.n)
            error('Incompatible sizes of vectors and matrix!');
        else
            input_type = 'Both are TTs';
        end
    end


    % Let v be the longer vector
    if strcmp(input_type, 'Both are vectors')
        if length(v) < length(w)
                [v, w] = deal(w, v);
        end
    end


    % Specify error calculation method
    if strcmp(input_type, 'Both are vectors')
        if ismatrix(err_type)
            method = 0;
        end
    elseif strcmp(input_type, 'Both are TTs')
        method = 0;
    end


    % Gauss quadrature rule
    if ~iscell(gqr)
        [pa, pw] = le_gqr(gqr(1));
        [qa, qw] = le_gqr(gqr(2));
        gqr = {[pa pw], [qa qw]};
    else
        [pa, pw] = deal(gqr{1}(:, 1), gqr{1}(:, 2));
        [qa, qw] = deal(gqr{2}(:, 1), gqr{2}(:, 2));
    end
    

    % Endpoints
    [xa, xb, ya, yb] = deal(xy_ep(1), xy_ep(2), xy_ep(3), xy_ep(4));

    
    % Case I: Both v and w are function handles
    if strcmp(input_type, 'Both are functions')
        % Normalization
        norm_sq_v = normsq_2d(v, xy_ep, err_type, gqr);
        if norm_sq_v ~= 0
            v = @(x, y)(v(x, y)/sqrt(norm_sq_v));
        end
        norm_sq_w = normsq_2d(w, xy_ep, err_type, gqr);
        if norm_sq_w ~= 0
            w = @(x, y)(w(x, y)/sqrt(norm_sq_w));
        end

        % Caculate the error
        [v, w] = deal(@(x, y)(v(x, y)+w(x, y)), @(x, y)(v(x, y)-w(x, y)));
        err = sqrt(min(normsq_2d(v, xy_ep, err_type, gqr), ...
                       normsq_2d(w, xy_ep, err_type, gqr)));
    end
    

    % Case II: v/w and w/v are a function handle and a vector resp.
    if strcmp(input_type, 'Function and vector')
        % Normalization
        norm_sq_v = normsq_2d(v, xy_ep, err_type, gqr);
        if norm_sq_v ~= 0
            v = @(x, y)(v(x, y)/sqrt(norm_sq_v));
        end
        norm_sq_w = normsq_2d(w, xy_ep, err_type, gqr);
        if norm_sq_w ~= 0
            w = w/sqrt(norm_sq_w);
        end

        % Pre-calculation
        N = sqrt(size(w, 1))+1; % Number of elements on each row
        M = N;  % Number of elements on each column
        hx = (xb-xa)/N;
        hy = (yb-ya)/M;
        J = hx*hy/4; % Jacobian
        w = reshape(w, [N-1 M-1])';
        w = [zeros(1, N+1); ...
             zeros(M-1, 1), w, zeros(M-1, 1); ...
             zeros(1, N+1)];

        % Shape function and gradients
        f = @(pn, qm, p, q)((1+pn*p)*(1+qm*q)/4);
        if strcmp(err_type, 'H1')
            fp = @(pn, qm, p, q)((pn*(1+qm*q))/4*2/hx);
            fq = @(pn, qm, p, q)(((1+pn*p)*qm)/4*2/hy);
            vx = matlabFunction(diff(sym(v), 'x'), 'Vars', {'x','y'});
            vy = matlabFunction(diff(sym(v), 'y'), 'Vars', {'x','y'});
        end

        % Calculate the error
        norm_sq_v = 0;
        norm_sq_w = 0;
        for m = 1:M % Each row of elements
        for n = 1:N % Each column of elements
            xn = xa+(n-1)*hx; % Left endpoint of x in the element
            ym = ya+(m-1)*hy; % Left endpoint of y in the element
            ww = @(p, q)(w(m, n)*f(-1, -1, p, q) ...
                         +w(m, n+1)*f(1, -1, p, q) ...
                         +w(m+1, n)*f(-1, 1, p, q) ...
                         +w(m+1, n+1)*f(1, 1, p, q));
            vv = @(p, q)(v((p+1)*hx/2+xn, (q+1)*hy/2+ym));
            if strcmp(err_type, 'H1')
                wp = @(p, q)(w(m, n)*fp(-1, -1, p, q) ...
                             +w(m, n+1)*fp(1, -1, p, q) ...
                             +w(m+1, n)*fp(-1, 1, p, q) ...
                             +w(m+1, n+1)*fp(1, 1, p, q));
                wq = @(p, q)(w(m, n)*fq(-1, -1, p, q) ...
                             +w(m, n+1)*fq(1, -1, p, q) ...
                             +w(m+1, n)*fq(-1, 1, p, q) ...
                             +w(m+1, n+1)*fq(1, 1, p, q));
                vp = @(p, q)(vx((p+1)*hx/2+xn, (q+1)*hy/2+ym));
                vq = @(p, q)(vy((p+1)*hx/2+xn, (q+1)*hy/2+ym));
                integrand_v = @(p, q)((ww(p, q)+vv(p, q))^2 ...
                                      +(wp(p, q)+vp(p, q))^2 ...
                                      +(wq(p, q)+vq(p, q))^2);
                integrand_w = @(p, q)((ww(p, q)-vv(p, q))^2 ...
                                      +(wp(p, q)-vp(p, q))^2 ...
                                      +(wq(p, q)-vq(p, q))^2);
            else
                integrand_v = @(p, q)((ww(p, q)+vv(p, q))^2);
                integrand_w = @(p, q)((ww(p, q)-vv(p, q))^2);
            end
            for k = 1:size(gqr{1}, 1)
            for l = 1:size(gqr{2}, 1)
                norm_sq_v = norm_sq_v+integrand_v(pa(k), qa(l))...
                                      *J*pw(k)*qw(l);
                norm_sq_w = norm_sq_w+integrand_w(pa(k), qa(l))...
                                      *J*pw(k)*qw(l);
            end
            end
        end
        end
        err = sqrt(min(norm_sq_v, norm_sq_w));
    end

    
    % Case III: Both v and w are vectors
    if strcmp(input_type, 'Both are vectors')
        % Interpolation
        v = reshape(v, [sqrt(numel(v)) sqrt(numel(v))])';
        w = reshape(w, [sqrt(numel(w)) sqrt(numel(w))])';
        v = [zeros(1, size(v, 2)+2); ...
             zeros(size(v, 1), 1), v, zeros(size(v, 1), 1); ...
             zeros(1, size(v, 2)+2)];
        w = [zeros(1, size(w, 2)+2); ...
             zeros(size(w, 1), 1), w, zeros(size(w, 1), 1); ...
             zeros(1, size(w, 2)+2)];
        [xv, yv] = meshgrid(linspace(xa, xb, size(v, 2)), ...
                            linspace(ya, yb, size(v, 1)));
        [xw, yw] = meshgrid(linspace(xa, xb, size(w, 2)), ...
                            linspace(ya, yb, size(w, 1)));
        w = interp2(xw, yw, w, xv, yv, 'linear');
        
        % Determine error type
        if method
            if ~ischar(err_type)
                method = 0;
            end
        else
            if ischar(err_type)
                if strcmp(err_type, 'L2')
                    err_type = fem_mat_2d_boost(size(v, 1)-1, xy_ep, ...
                                                'Mass', gqr);
                elseif strcmp(err_type, 'H1')
                    err_type = fem_mat_2d_boost(size(v, 1)-1, xy_ep, ...
                                                'Laplace+Mass', gqr);
                else
                    error(['If norm_type is a character, ', ...
                           'it should be L2/H1']);
                end
            end
        end

        % Normalization
        norm_sq_v = normsq_2d(v, xy_ep, err_type, gqr, method);
        if norm_sq_v ~= 0
            v = v/sqrt(norm_sq_v);
        end
        norm_sq_w = normsq_2d(w, xy_ep, err_type, gqr, method);
        if norm_sq_w ~= 0
            w = w/sqrt(norm_sq_w);
        end

        % Calculate the error
        [v, w] = deal(v+w, v-w);
        norm_sq_v = normsq_2d(v, xy_ep, err_type, gqr, method);
        norm_sq_w = normsq_2d(w, xy_ep, err_type, gqr, method);
        err = sqrt(min(norm_sq_v, norm_sq_w));
    end

    
    % Case IV: Both v and w are TT-tensors
    if strcmp(input_type, 'Both are TTs')
        [v, w] = deal(round(v+w, tol), round(v-w, tol));
        if ischar(err_type)
            if strcmp(err_type, 'L2')
                err_type = tt_qfemmat2d(v.d/2, ...
                                        xy_ep([2 4])-xy_ep([1 3]), 'Mass');
            else
                err_type = tt_qfemmat2d(v.d/2, ...
                                        xy_ep([2 4])-xy_ep([1 3]), ...
                                        'Laplace+Mass');
                
            end
        end
        err = min(dot(v, round(err_type*v, tol)), ...
                  dot(w, round(err_type*w, tol)));
    end


end