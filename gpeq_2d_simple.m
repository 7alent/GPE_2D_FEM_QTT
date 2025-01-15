% 
% Solve 2-D Gross–Pitaevskii Equation via SCF Iteration
% 
% [U, D] = GPEQ_2D_SIMPLE(N, VARARGIN)
%   Solve 2-D Gross–Pitaevskii Equation:
% 
%                   [-(1/2)△+v(x)+u(x)^2]u(x) = λu(x)
%   
%   with Dirichlet boundary on 2-D square area using FEM-SCF iteration
% 
%   [Input Argument]
%       N - Integer, number of elements on X/Y axis, total number of
%           elements will be N^2
%       xy_ep - Vector, optional, set it as [xa xb ya yb] if the domain is 
%               [xa, xb] x [ya, yb] (default: [0 1 0 1])
%       v - Function handle or sym, optional, the potential function 
%           (default: [])
%       gqr - Vector or cell, optional, set it as a vector [po qo] to set 
%             the orders of Legendre Gauss quadrature of local integral as 
%             po/qo corresponding to local variables p/q (w.r.t. global 
%             variables x/y), or a 2 x 1 cell whose entries are 2-column 
%             matrices (the 1-st column is the abscissas while the 2-nd is 
%             the weights) containing the quadrature rules of X/Y axis 
%             resp. (default: [5 5])
%       tol - Scalar, optional, tolerance of the iteration and eigen solver
%             (default: 1e-14)
%       max_iter - Scalar, optional, maximum iteration times of SCF 
%                  iteration and eigen solver (default: 100)
%       U0 - Vector, optional, initial guess of eigen function (default: 
%            random numbers from standard normal distribution)
%       method - Cell or character, optional, which method to use for
%                solving the eigen problem, should be a subset of {'FEM', 
%                'ALS', 'DMRG'} (eigs(), 1-site DMRG and 2-site DMRG resp., 
%                the last 2 methods are QTT methods and require N = 2^d+1)
%                (default: {'FEM'})
%       rmax - Scalar, optional, max TT rank (default: 8)
% 
%   [Ouput Argument]
%       output - Structure with fields:
%           eigen_function - Cell, eigen vectors
%           eigen_value - Table, eigen values
%           H1_mat - Cell, sum of Laplace and mass matrices (and their 
%                    QTT), this can be used to calculate H1 norm
%           iter_info - Table, more details of iteration, contains 
%                        rows:
%               copnverge - Boolen, whether the iteration converges 
%                           (1: Yes, 0: No)
%               residual of rho - Scalar, residual (e.g. H1-error) of 
%                                 electron density (function)
%               total time - Scalar, total time of SCF iteration (unit: 
%                            second)
%               iteration - Scalar, total iterations
% 
% Details:
%   1. In output, results will be in the order as {'FEM', 'ALS', 'DMRG'}
%      rather than the order of method


function output = gpeq_2d_simple(N, varargin)
    %% Input assignment and type check

    %%% Intialization
    [xy_ep, v, gqr, tol, max_iter, ... 
     U0, U0_als, U0_dmrg, method, rmax] = deal([]);


    %%% Assignment
    for i = 1:(nargin-1)

        % Numeric inputs
        if isvector(varargin{i}) && isnumeric(varargin{i})
            if length(varargin{i}) == 4
                xy_ep = varargin{i};
            elseif length(varargin{i}) == 2
                gqr = varargin{i};
            elseif isscalar(varargin{i}) % Numeric inputs will be assigned 
                                         % to max_iter, rmax and tol 
                                         % successively
                if rem(varargin{i}, 1) == 0
                    if isempty(max_iter)
                        max_iter = varargin{i};
                    else
                        rmax = varargin{i};
                    end
                else
                    tol = varargin{i};
                end
            else
                U0 = varargin{i};
            end

        % Function input
        elseif isa(varargin{i}, 'function_handle')
            v = varargin{i};
        elseif isa(varargin{i}, 'sym')
            v = matlabFunction(varargin{i});

        % Character and cell inputs
        elseif ischar(varargin{i})
            method = cell(1);
            method{1} = varargin{i};
        elseif iscell(varargin{i})
            method = varargin{i};
        else
            error('Wrong input!');
        end

    end


    %%% Check number of element
    if ~(isnumeric(N) && isscalar(N))
        error('Number of element on X/Y axis should be a scalar!');
    elseif rem(N, 1) ~= 0 || N <= 0
        error(['Number of element on X/Y axis should be a positive ', ...
               'integer!']);
    end
    

    %%% Check methods
    if ~isempty(method)
        if ~any(strcmp(method, 'FEM')) && ~any(strcmp(method, 'ALS')) ... 
           && ~any(strcmp(method, 'DMRG'))
            error(['Methods should be a subset of {''FEM'', ''ALS'', ', ...
                   '''DMEG''}!']);
        end
    end


    %%% Sort methods
    if ~isempty(method)
        method = char_sort(method, {'FEM', 'ALS', 'DMRG'});
    end


    %%% Default input
    if isempty(xy_ep)
        xy_ep = [0 1 0 1];
    end
    if isempty(gqr)
        gqr = [5 5];
    end
    if isempty(tol)
        tol = 1e-14;
    end
    if isempty(max_iter)
        max_iter = 100;
    end
    if isempty(method)
        method = {'FEM'};
    end

    
    %%% Calculate level of the mesh
    if any(strcmp(method, 'ALS')) || any(strcmp(method, 'DMRG'))
        d = log(N-1)/log(2);
        if rem(d, 1) ~= 0
            error('Level of the mesh should be an integer!');
        end
    end


    %%% Calculate domain sizes and spacing
    dom_size = xy_ep([2 4])-xy_ep([1 3]);
    h_x = dom_size(1)/N;
    h_y = dom_size(2)/N;


    %% Generate start vectors

    if isempty(U0)
        % Default value for max TT rank
        if isempty(rmax)
            rmax = 8;
        end
        
        % Start vector of ALS
        if any(strcmp(method, 'ALS'))
            U0_als = tt_rand(2, 2*d, rmax);
        end
        
        % Start vector of 2-site DMRG
        if any(strcmp(method, 'DMRG'))
            if any(strcmp(method, 'ALS')) % We can use start vector of ALS
                U0_dmrg = U0_als;
            else
                U0_dmrg = tt_rand(2, 2*d, rmax);
            end
        end
        
        % Start vector of FEM
        if any(strcmp(method, 'FEM'))
            if any(strcmp(method, 'ALS')) % We can use start vector of ALS
                U0 = full(U0_als);
            elseif any(strcmp(method, 'DMRG')) % Same as above
                U0 = full(U0_dmrg);
            else
                U0 = randn((N-1)^2, 1);
            end
        end

    else

        % U0 is a TT-vector
        if ~isvector(U0)
            if any(strcmp(method, 'ALS'))
                U0_als = U0;
            end
            if any(strcmp(method, 'DMRG'))
                U0_dmrg = U0;
            end
            if any(strcmp(method, 'FEM'))
                U0 = full(U0);
            end

        % U0 is a vector
        else
            if any(strcmp(method, 'ALS'))
                U0_als = tt_tensor(tt_qfromfull(U0, 1, 2*d, tol, 2));
            end
            if any(strcmp(method, 'DMRG'))
                if any(strcmp(method, 'ALS'))
                    U0_dmrg = U0_als; % We can use start vector of ALS
                else
                    U0_dmrg = tt_tensor(tt_qfromfull(U0, 1, 2*d, tol, 2));
                end
            end
        end

    end

    
    %% Assemle stiffness matrices

    %%% Assemble Laplace and mass stiffness matrices
    fprintf('\nAssembling Laplace and mass stiffness matrices... ');
    tic;
    [A, M, K] = fem_mat_2d_exact(N, xy_ep([2 4])-xy_ep([1 3]));
    A = A/2;
    t = toc;
    fprintf('finished (Time: %.4fs)\n', t);
    

    %%% Assemble QTT of Laplace and mass stiffness matrices
    if any(strcmp(method, 'ALS')) || any(strcmp(method, 'DMRG'))
        fprintf(['\nCalculating QTT of Laplace and mass stiffness ', ...
                 'matrices... ']);
        tic;
        [A_qtt, M_qtt, K_qtt] = tt_qfemmat2d(d, xy_ep([2 4])-xy_ep([1 3]));
        A_qtt = A_qtt/2;
        t = toc;
        fprintf('finished (Time: %.4fs)\n', t);
    end
    
    %%% Assemble potential stiffness matrix
    if ~isempty(v)
        fprintf('Assembling external potential stiffness matrix... ');
        tic;
        B = fem_mat_2d_boost(v, N, xy_ep, 'Potential', ...
                             gqr); % External potential stiffness matrix
        t = toc;
        fprintf('finished (Time: %.4fs)\n', t);

        % Assemble QTT of potential stiffness matrix
        if any(strcmp(method, 'ALS')) || any(strcmp(method, 'DMRG'))
            fprintf(['\nCalculating QTT of external potential ', ...
                     'stiffness matrices... ']);
            tic;
            B_qtt = tt_qband(B, tol);
            t = toc;
            fprintf('finished (Time: %.4fs)\n', t);
        end

    else
        fprintf(['No external potential given, skip potential ', ...
                 'stiffness matrix calculation\n']);
    end


    %% Generate start values (FEM)

    if any(strcmp(method, 'FEM'))
        fprintf('\nStart initialization for FEM...\n');
        
        % Assemble stiffness matrix for electron density
        fprintf('Assembling stiffness matrix of u(x)^2... ');
        tic;
        C = sepdiag2d(N-1, [36 2 6 1]*h_x*h_y/144, 1)*(U0.^2); % Stiffness 
                                                               % matrix of 
                                                               % u(x)^2
        C = reshape(C, [(N-1)^2 (N-1)^2]);
        t = toc;
        fprintf('finished (Time: %.4fs)\n', t);
        
        % Solve eigen problem
        fprintf('Solving eigen problem... ');
        tic;
        if ~isempty(v)
            L = A+B+C;
        else
            L = A+C;
        end
        [U, D] = eigs(L, M, 1, 'smallestabs', 'Tolerance', tol, ...
                      'StartVector', U0, 'MaxIterations', max_iter);
        rho_res = fem_err_2d_tt(U.^2, U0.^2, K); % Residual of u(x)^2
        t = toc;
        fprintf('finished (Time: %.4fs)\nInitialization completed!\n', t);

    end
    
    
    %% Generate start values (ALS)

    if any(strcmp(method, 'ALS'))
        fprintf('\nStart initialization for ALS...\n');
        
        % Assemble stiffness matrix for electron density
        fprintf('Assembling stiffness matrix of u(x)^2... ');
        tic;
        C_als = tt_qsolpot2d(round(U0_als.^2, tol), xy_ep, ...
                             tol); % Stiffness matrix of u(x)^2
        t = toc;
        fprintf('finished (Time: %.4fs)\n', t);
        
        % Solve eigen problem
        fprintf('Solving eigen problem... ');
        tic;
        if ~isempty(v)
            L_als = round(round(A_qtt+B_qtt, tol)+C_als, tol);
        else
            L_als = round(A_qtt+C_als, tol);
        end
        [U_als, D_als, ~] = dmrg2(L_als, M_qtt, tol, 'b', 1, ...
                                  'numblocks', 1, ...
                                  'local_iters', max_iter, ...
                                  'rmax', rmax, 'verb', 0, 'x0', U0_als);
        rho_res_als = fem_err_2d_tt(round(U_als.^2, tol), ...
                                    round(U0_als.^2, tol), ...
                                    K_qtt, tol); % Residual of u(x)^2
        t = toc;
        fprintf('finished (Time: %.4fs)\nInitialization completed!\n', t);

    end
    

    %% Generate start values (DMRG)

    if any(strcmp(method, 'DMRG'))
        fprintf('\nStart initialization for 2-site DMRG...\n');
        
        % Assemble stiffness matrix for electron density
        fprintf('Assembling stiffness matrix of u(x)^2... ');
        tic;
        C_dmrg = tt_qsolpot2d(round(U0_dmrg.^2, tol), xy_ep, d, ...
                              tol); % Stiffness matrix of u(x)^2
        t = toc;
        fprintf('finished (Time: %.4fs)\n', t);
        
        % Solve eigen problem
        fprintf('Solving eigen problem... ');
        tic;
        if ~isempty(v)
            L_dmrg = round(round(A_qtt+B_qtt, tol)+C_dmrg, tol);
        else
            L_dmrg = round(A_qtt+C_dmrg, tol);
        end
        [U_dmrg, D_dmrg, ~] = dmrg2(L_dmrg, M_qtt, tol, 'b', 1, ...
                                    'numblocks', 2, ...
                                    'local_iters', max_iter, ...
                                    'rmax', rmax, 'verb', 0, ...
                                    'x0', U0_dmrg);
        rho_res_dmrg = fem_err_2d_tt(round(U_dmrg.^2, tol), ...
                                     round(U0_dmrg.^2, tol), ...
                                     K_qtt, tol); % Residual of u(x)^2
        t = toc;
        fprintf('finished (Time: %.4fs)\nInitialization completed!\n', t);

    end


    %% SCF Iteration (FEM)

    if any(strcmp(method, 'FEM'))

        %%% Intialization
        iter = 0;
        total_time = 0;
        fprintf('\nStart SCF iteration of FEM...\n');
        
        %%% SCF
        while rho_res >= 10*tol && iter < max_iter
    
            % Iteration counter
            iter = iter+1;
    
            tic;
    
            % Mix electron density
            U0 = U;
            C = sepdiag2d(N-1, [36 2 6 1]*h_x*h_y/144, ...
                          1)*(U0.^2); % Stiffness matrix of u(x)^2
            C = reshape(C, [(N-1)^2 (N-1)^2]);
            
            % Calculate left-hand-side matrix
            if ~isempty(v)
                L = A+B+C;
            else
                L = A+C;
            end
    
            % Solve eigen problem
            [U, D] = eigs(L, M, 1, 'smallestabs', 'Tolerance', tol, ...
                          'StartVector', U0, 'MaxIterations', max_iter);
    
            % Calculate error of electron density
            rho_res = fem_err_2d_tt(U.^2, U0.^2, K); % Residual of u(x)^2
    
            t = toc;
    
            fprintf(['\nIteration: %d \t Eigen Value: %.8f \t ', ...
                     'Residual: %.4e \t Time: %.4fs'], ...
                    iter, D, rho_res, t);
    
            % Iteration timer
            total_time = total_time+t;
        end
    
        
        %%% Convergence check (FEM)
        if rho_res < 10*tol
            converge = 1;
            fprintf(['\n\nFEM converges successfully ', ...
                     '(Total time: %.4fs)!\n\n'], total_time);
        else
            converge = 0;
            fprintf(['\n\nFEM fails to converge (Total time: %.4fs)!', ...
                     '\n\n'], total_time);
        end


    end

    
    %% SCF Iteration (ALS)

    if any(strcmp(method, 'ALS'))

        %%% Intialization
        iter_als = 0;
        total_time_als = 0;
        fprintf('\nStart SCF iteration of ALS...\n');
        
        %%% SCF
        while rho_res_als >= 10*tol && iter_als < max_iter
            
            % Iteration counter
            iter_als = iter_als+1;
    
            tic;
    
            % Mix electron density
            U0_als = U_als;
            C_als = tt_qsolpot2d(round(U0_als.^2, tol), xy_ep, d, ...
                                 tol); % Stiffness matrix of u(x)^2
    
            % Calculate left-hand-side matrix
            if ~isempty(v)
                L_als = round(round(A_qtt+B_qtt, tol)+C_als, tol);
            else
                L_als = round(A_qtt+C_als, tol);
            end
            
            % Solve eigen problem
            [U_als, D_als, ~] = dmrg2(L_als, M_qtt, tol, 'b', 1, ...
                                      'numblocks', 1, ...
                                      'local_iters', max_iter, ...
                                      'rmax', rmax, 'verb', 0, ...
                                      'x0', U0_als);
    
            % Calculate error of electron density
            rho_res_als = fem_err_2d_tt(round(U_als.^2, tol), ...
                                        round(U0_als.^2, tol), ...
                                        K_qtt, tol); % Residual of u(x)^2
            
            t = toc;
    
            fprintf(['\nIteration: %d \t Eigen Value: %.8f \t ', ...
                     'Residual: %.4e \t Time: %.4fs'], ...
                    iter_als, D_als, rho_res_als, t);
    
            % Iteration timer
            total_time_als = total_time_als+t;
        end
    
        
        %%% Convergence check (ALS)
        if rho_res_als < 10*tol
            converge_als = 1;
            fprintf(['\n\nALS converges successfully ', ...
                     '(Total time: %.4fs)!\n\n'], total_time_als);
        else
            converge_als = 0;
            fprintf(['\n\nALS fails to converge (Total time: %.4fs)!', ...
                     '\n\n'], total_time_als);
        end


    end


    %% SCF Iteration (2-site DMRG)

    if any(strcmp(method, 'DMRG'))

        %%% Intialization
        iter_dmrg = 0;
        total_time_dmrg = 0;
        fprintf('\nStart SCF iteration of 2-site DMRG...\n');
        
        %%% SCF
        while rho_res_dmrg >= 10*tol && iter_dmrg < max_iter
    
            % Iteration counter
            iter_dmrg = iter_dmrg+1;
    
            tic;
    
            % Mix electron density
            U0_dmrg = U_dmrg;
            C_dmrg = tt_qsolpot2d(round(U0_dmrg.^2, tol), xy_ep, d, ...
                                  tol); % Stiffness matrix of u(x)^2
    
            % Calculate left-hand-side matrix
            if ~isempty(v)
                L_dmrg = round(round(A_qtt+B_qtt, tol)+C_dmrg, tol);
            else
                L_dmrg = round(A_qtt+C_dmrg, tol);
            end
    
            % Solve eigen problem
            [U_dmrg, D_dmrg, ~] = dmrg2(L_dmrg, M_qtt, tol, 'b', 1, ...
                                        'numblocks', 2, ...
                                        'local_iters', max_iter, ...
                                        'rmax', rmax, 'verb', 0, ...
                                        'x0', U0_dmrg);
    
            % Calculate error of electron density
            rho_res_dmrg = fem_err_2d_tt(round(U_dmrg.^2, tol), ...
                                         round(U0_dmrg.^2, tol), ...
                                         K_qtt, tol); % Residual of u(x)^2
            
            t = toc;
            
            fprintf(['\nIteration: %d \t Eigen Value: %.8f \t ', ...
                     'Residual: %.4e \t Time: %.4fs'], ...
                    iter_dmrg, D_dmrg, rho_res_dmrg, t);
    
            % Iteration timer
            total_time_dmrg = total_time_dmrg+t;
    
        end
    
        
        %%% Convergence check (DMRG)
        if rho_res_dmrg < 10*tol
            converge_dmrg = 1;
            fprintf(['\n\n2-site DMRG converges successfully ', ...
                     '(Total time: %.4fs)!\n\n'], total_time_dmrg);
        else
            converge_dmrg = 0;
            fprintf(['\n\n2-site DMRG fails to converge ', ...
                     '(Total time: %.4fs)!\n\n'], total_time_dmrg);
        end


    end


    %% Return output

    %%% Initialize outputs
    U_struct = struct();
    D_table = table('Size', [1 length(method)], 'VariableTypes', ...
                    repmat({'double'}, [length(method) 1]), ...
                    'VariableNames', method);
    info_table = table('Size', [4 length(method)], 'VariableTypes', ...
                       repmat({'double'}, [length(method) 1]), ...
                       'VariableNames', method, 'RowNames', ...
                       {'converge', 'residual of rho', 'total time', ...
                        'iteration'});


    %%% Pack the sum of Laplace and mass matrices
    if any(strcmp(method, 'ALS')) || any(strcmp(method, 'DMRG'))
        K_cell = {K, K_qtt};
    else
        K_cell = K;
    end
    

    %%% Pack solution and iteration information
    for i = 1:length(method)
        switch method{i}
            case 'FEM'
                U_struct.FEM = U;
                D_table(1, i) = {D};
                info_table('converge', i) = {converge};
                info_table('residual of rho', i) = {rho_res};
                info_table('total time', i) = {total_time};
                info_table('iteration', i) = {iter};
            case 'ALS'
                U_struct.ALS = U_als;
                D_table(1, i) = {D_als};
                info_table('converge', i) = {converge_als};
                info_table('residual of rho', i) = {rho_res_als};
                info_table('total time', i) = {total_time_als};
                info_table('iteration', i) = {iter_als};
            case 'DMRG'
                U_struct.DMRG = U_dmrg;
                D_table(1, i) = {D_dmrg};
                info_table('converge', i) = {converge_dmrg};
                info_table('residual of rho', i) = {rho_res_dmrg};
                info_table('total time', i) = {total_time_dmrg};
                info_table('iteration', i) = {iter_dmrg};
        end
    end


    %%% Pack the final output
    output = struct();
    output.eigen_function = U_struct;
    output.eigen_value = D_table;
    output.H1_mat = K_cell;
    output.iter_info = info_table;


end