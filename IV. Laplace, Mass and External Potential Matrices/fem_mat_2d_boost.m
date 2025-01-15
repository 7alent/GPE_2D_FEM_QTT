% 
% 2-D Finite Element Stiffness Matrices
% 
% VARARGOUT = FEM_MAT_2D_BOOST(VARARGIN)
%   Assemble (Laplace, potential and mass) stiffness matrices in 2-D finite 
%   element method
%   The algorithm will be boosted by assembling upper triangular parts of 
%   the matrices and using symmetry to get the whole matrices
% 
%   [Input Argument]
%       v - Function handle or sym, optional, the potential function
%       N - Integer, number of elements on X/Y axis, total number of
%           elements will be N^2
%       xy_bound - Vector, set it as [xa xb ya yb] if the domain is 
%                  [xa, xb] x [ya, yb]
%       mat_type - Cell, a subset of {'Laplace', 'Potential', 'Mass', 
%                  'Laplace+Potential', 'Laplace+Mass'} that denotes types 
%                  of stiffness matrices to be assembled
%       gqr - Vector or cell, set it as a vector [po qo] to set the orders 
%             of Legendre Gauss quadrature of local integral as po/qo 
%             corresponding to local variables p/q (w.r.t. global variables 
%             x/y), or a 2 x 1 cell whose entries are 2-column matrices 
%             (the 1-st column is the abscissas while the 2-nd is the
%             weights) containing the quadrature rules of X/Y axis resp.
% 
%   [Ouput Argument]
%       A, B, M, H, K - Matrix, the first 3 outputs denotes Laplace, 
%                       potential, mass matrices resp., while the last 2
%                       outputs are (A+B) and (A+M) in fact, outputs will
%                       be in the order as mat_type
% 
% Details:
%   1. The stiffness matrices originate from FEM for solving 2-D elliptical
%      PDE eigen problem:
% 
%                           -Δu + v * u = λ * u
% 
%      with u = u(x, y) is the solution (eigen function) and v = v(x, y) is
%      the potential function, where (x, y) ∈ [xa, xb] x [ya, yb]
%   2. Some abbreviations for object names:
%                   N - Node
%                   E - Element
%                   n - (total) number
%                   gi - global index
%                   li - local index
%                   mi - (stiffness) matrix (entry) index
%                   Tr - trial function
%                   Te - test function
%                   gc/lc - global/local coordinate
%                   ip - interior point(not a boundary point)
%                   a - (Gauss point) abscissa 
%                   w - (Gauss) weight
%   3. All nodes (locally/globally) are indexed from left to the right and
%      from bottom to the top; indices of entries of stiffness matrices are
%      in fact the indices of interior nodes (points) and they are indexed
%      similarly; a simple examle is the case N = 3:
%                           21  22  23  24  25
%                           16  17  18  19  20
%                           11  12  13  14  15
%                           6   7   8   9   10
%                           1   2   3   4   5
%      global indices are shown above; when we focus on an element (the 
%      first element, for example), local indices of 4 nodes of it (i.e. 
%      the 1, 2, 6, 7-th nodes) are always 1~4; interior points are 7~9, 
%      12~14 and 17~19-th nodes and their matrix indices are 1~9
%   4. Both trial and test functions are piecewise linear functions
%   5. The output matrix H is recommended to use in the solution of eigen 
%      problem 
% 
%                              H * u = λ * M * u
% 
%      derived from the original problem (note that here u is the weight 
%      vector of basis functions);
%      the output matrix M/K is used to caculate the L2/H1 norm of u:
% 
%                       norm_{L2}(u) = sqrt(u' * M * u)
%                       norm_{H1}(u) = sqrt(u' * K * u)
% 
%      it's faster to caculate H or K directly than to caculate (A and B) 
%      or (A and M) separately and add them up
%   6. When the potential function has singular points, e.g. 
% 
%                           v(x, y) = 1/abs(x-y)
% 
%      orders of Gauss quadrature should be chosen so that Gauss points
%      don't fall into any singular point, e.g. if we choose [5 5] then
%      Gauss points will be on the singular line y = x, so we'd better use
%      different orders for x/y like [5 4]
%   7. To generate A, M, H and K, codes below are faster:
% 
%           dom_size = xy_bound([2 4])-xy_bound([1 3]);
%           A = fem_mat_2d_exact(N, dom_size, 'Laplace');
%           M = fem_mat_2d_exact(N, dom_size, 'Mass');
%           H = A+fem_mat_2d_boost(v, N, xy_bound, [po qo]);
%           K = fem_mat_2d_exact(N, dom_size, 'Laplace+Mass');


function varargout = fem_mat_2d_boost(varargin)
    % Input assignment
    if nargin == 4
        v = [];
        N = varargin{1};
        xy_bound = varargin{2};
        mat_type = varargin{3};
        gqr = varargin{4};
        if any(strcmp(mat_type, 'Potential')) || ...
           any(strcmp(mat_type, 'Laplace+Potential'))
            error(['When calculating potential stiffness matrix, ', ...
                   'potential function should be given!']);
        end
    elseif nargin == 5
        v = varargin{1};
        N = varargin{2};
        xy_bound = varargin{3};
        mat_type = varargin{4};
        gqr = varargin{5};
        if (any(strcmp(mat_type, 'Potential')) || ...
            any(strcmp(mat_type, 'Laplace+Potential'))) && isempty(v)
            error(['When calculating potential stiffness matrix, ', ...
                   'potential function should be given!']);
        end
    else
        error('Number of arguments should be 4 or 5!');
    end


    % Input check
    if ~iscell(mat_type)
        mat_type = {mat_type};
    end
    if ~isempty(v) && ~isa(v, 'function_handle')
        v = matlabFunction(v);
    end
    

    % Step length
    hx = (xy_bound(2)-xy_bound(1))/N;
    hy = (xy_bound(4)-xy_bound(3))/N;
    

    % Global coordinates of nodes
    N_n = (N+1)^2; % Number of nodes
    x0 = xy_bound(1);
    y0 = xy_bound(3);
    N_gc = zeros(2, N_n); % each column corresponds to a node
    for i = 1:(N+1) % i-th row of nodes
        for j = 1:(N+1) % j-th column of nodes
            N_gc(:, (i-1)*(N+1)+j) = [x0+(j-1)*hx; y0+(i-1)*hy];
        end
    end
    

    % Global indices of nodes in each element
    % Type (interior/boundary) of nodes
    E_n = N^2; % Number of elements
    EN_gi = zeros(4, E_n); % Its n-th column contains global indices of all
                           % nodes in the n-th element
    ip = ones(N_n, 1); % Whether the n-th node is an interior point, 
                       % if it is: set the n-th entry as its matrix index; 
                       % if not: set it as 0
    for i = 1:N % i-th row of elements
        for j = 1:N % j-th column of elements
            left_bottom_node = (i-1)*(N+1)+j;
            EN_gi(:, (i-1)*N+j) = [left_bottom_node; ...
                                   left_bottom_node+1; ...
                                   left_bottom_node+N+1; ...
                                   left_bottom_node+N+2];
            if i == 1 % Deal with boundary nodes
                ip(EN_gi([1 2], (i-1)*N+j)) = 0;
            elseif i == N
                ip(EN_gi([3 4], (i-1)*N+j)) = 0;
            end
            if j == 1
                ip(EN_gi([1 3], (i-1)*N+j)) = 0;
            elseif j == N
                ip(EN_gi([2 4], (i-1)*N+j)) = 0;
            end
        end
    end
    mi = 0; % (Number of) matrix indices
    for i = 1:N_n
        if ip(i)
            mi = mi+1;
            ip(i) = mi;
        end
    end
    

    % Nodes and number of (local basis) trial/test functions
    ETrN_gi = EN_gi; % Its n-th column contains global indices of all nodes 
                     % w.r.t. each trial function in the n-th element
    ETeN_gi = EN_gi; % Its n-th column contains global indices of all nodes 
                     % w.r.t. each test function in the n-th element
    Tr_n = size(ETrN_gi, 1); % Number of trial functions
    Te_n = size(ETeN_gi, 1); % Number of test functions
    

    % Intiallize output
    if any(strcmp(mat_type, 'Laplace')) % Laplace matrix
        A = sparse(mi, mi);
    end
    if any(strcmp(mat_type, 'Potential')) % Potential matrix
        B = sparse(mi, mi);
    end
    if any(strcmp(mat_type, 'Mass')) % Mass matrix
        M = sparse(mi, mi);
    end
    if any(strcmp(mat_type, 'Laplace+Potential')) % Sum of Laplace and 
                                                  % potential matrix
        H = sparse(mi, mi);
    end
    if any(strcmp(mat_type, 'Laplace+Mass')) % Sum of Laplace and mass 
                                             % matrix
        K = sparse(mi, mi);
    end
    

    % Local basis and their gradients
    syms p q % Local variables
    syms p_i q_i % Coordinate of local nodes
    syms f(p_i, q_i, p, q) % (Piecewise linear) shape function
    f(p_i, q_i, p, q) = (1+p_i*p)*(1+q_i*q)/4;
    if any(strcmp(mat_type, 'Laplace')) || ...
       any(strcmp(mat_type, 'Laplace+Potential')) || ...
       any(strcmp(mat_type, 'Laplace+Mass'))
        fp = matlabFunction(diff(f, p)); % Gradients of shape function
        fq = matlabFunction(diff(f, q));
    end
    f = matlabFunction(f);
    EN_lc = [[-1; -1], [1; -1], [-1; 1], [1; 1]]; % Local coordinates of 
                                                  % all nodes of an element
    
    
    % Jacobian in change of variables formula of integral
    J = hx*hy/4;


    % Gauss quadrature rules
    if iscell(gqr)
        [pa, pw] = deal(gqr{1}(:, 1), gqr{1}(:, 2));
        [qa, qw] = deal(gqr{2}(:, 1), gqr{2}(:, 2));
    else
        [pa, pw] = le_gqr(gqr(1));
        [qa, qw] = le_gqr(gqr(2));
    end


    % Assemble stiffness matrices
    for n = 1:E_n % Focus on each element at each step

        % Mind the multiples of basis functions and their gradients (which 
        % originate from change of variables formula) !
        for i = 1:Tr_n % Number of (nodes of) trial functions
        for j = 1:Te_n % Number of (nodes of) test functions

            % Only interior points are considered
            % Only calculate upper triangular parts of matrices
            if ip(ETrN_gi(i, n)) && ip(ETeN_gi(j, n)) && ...
               ip(ETrN_gi(i, n)) <= ip(ETeN_gi(j, n))
                % Trial/test function and their gradients
                Tr = @(p, q)(f(EN_lc(1, i), EN_lc(2, i), p, q));
                Te = @(p, q)(f(EN_lc(1, j), EN_lc(2, j), p, q));
                if any(strcmp(mat_type, 'Laplace')) || ...
                   any(strcmp(mat_type, 'Laplace+Potential')) || ...
                   any(strcmp(mat_type, 'Laplace+Mass'))
                    % Multiples (2/hx and 2/hy) of trial/test functions 
                    % originate from chain rule of differentiation: 
                    %               ∂g/∂x = (∂g/∂p) * (dp/dx)
                    %               ∂g/∂y = (∂g/∂q) * (dq/dy)
                    Trp = @(p, q)(fp(EN_lc(1, i), EN_lc(2, i), p, q)*2/hx);
                    Trq = @(p, q)(fq(EN_lc(1, i), EN_lc(2, i), p, q)*2/hy);
                    Tep = @(p, q)(fp(EN_lc(1, j), EN_lc(2, j), p, q)*2/hx);
                    Teq = @(p, q)(fq(EN_lc(1, j), EN_lc(2, j), p, q)*2/hy);
                end
                
                % Assemble Laplacian
                if any(strcmp(mat_type, 'Laplace'))
                    int_A = @(p, q)...
                            ((Trp(p, q)*Tep(p, q)+Trq(p, q)*Teq(p, q))*J);
                    for k = 1:length(pa)
                    for l = 1:length(qa)
                        A(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) = ...
                        A(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) + ...
                        int_A(pa(k), qa(l))*pw(k)*qw(l);
                    end
                    end
                end
                
                % Assemble potential matrix
                if any(strcmp(mat_type, 'Potential'))
                    for k = 1:length(pa)
                    for l = 1:length(qa)
                        s = v((pa(k)+1)*hx/2+N_gc(1, EN_gi(1, n)), ...
                              (qa(l)+1)*hy/2+N_gc(2, EN_gi(1, n)));
                        int_B = @(p, q)(s*Tr(p, q)*Te(p, q)*J);
                        B(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) = ...
                        B(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) + ...
                        int_B(pa(k), qa(l))*pw(k)*qw(l);
                    end
                    end
                end
                
                % Assemble mass matrix
                if any(strcmp(mat_type, 'Mass'))
                    int_M = @(p, q)(Tr(p, q)*Te(p, q)*J);
                    for k = 1:length(pa)
                    for l = 1:length(qa)
                        M(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) = ...
                        M(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) + ...
                        int_M(pa(k), qa(l))*pw(k)*qw(l);
                    end
                    end
                end

                % Assemble the sum of Laplacian and potential matrix
                if any(strcmp(mat_type, 'Laplace+Potential'))
                    int_A = @(p, q)...
                            ((Trp(p, q)*Tep(p, q)+Trq(p, q)*Teq(p, q))*J);
                    for k = 1:length(pa)
                    for l = 1:length(qa)
                        s = v((pa(k)+1)*hx/2+N_gc(1, EN_gi(1, n)), ...
                              (qa(l)+1)*hy/2+N_gc(2, EN_gi(1, n)));
                        int_B = @(p, q)(s*Tr(p, q)*Te(p, q)*J);
                        int_H = @(p, q)(int_A(p, q)+int_B(p, q));
                        H(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) = ...
                        H(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) + ...
                        int_H(pa(k), qa(l))*pw(k)*qw(l);
                    end
                    end
                end

                % Assemble the sum of Laplacian and mass matrix
                if any(strcmp(mat_type, 'Laplace+Mass'))
                    int_A = @(p, q)...
                            ((Trp(p, q)*Tep(p, q)+Trq(p, q)*Teq(p, q))*J);
                    int_M = @(p, q)(Tr(p, q)*Te(p, q)*J);
                    int_K = @(p, q)(int_A(p, q)+int_M(p, q));
                    for k = 1:length(pa)
                    for l = 1:length(qa)
                        K(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) = ...
                        K(ip(ETrN_gi(i, n)), ip(ETeN_gi(j, n))) + ...
                        int_K(pa(k), qa(l))*pw(k)*qw(l);
                    end
                    end
                end
                    
            end

        end
        end

    end


    % Output
    varargout = cell(length(mat_type), 1);
    for i = 1:length(mat_type)
        switch mat_type{i}
            case 'Laplace'
                varargout{i} = A+triu(A, 1)';
            case 'Potential'
                varargout{i} = B+triu(B, 1)';
            case 'Mass'
                varargout{i} = M+triu(M, 1)';
            case 'Laplace+Potential'
                varargout{i} = H+triu(H, 1)';
            case 'Laplace+Mass'
                varargout{i} = K+triu(K, 1)';
        end
    end
    
    
end