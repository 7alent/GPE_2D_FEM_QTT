% 
% Calculate Legendre Gauss Quadrature
% 
% output = GQ2(INTEGRAND, INTERVAL, RULE)
%   Given a 2-D function, calulate its double integral on a certain area
%   via Legendre Gauss Quadrature
% 
%   [Input Argument]
%       integrand - Function handle or sym, the function to be integrated
%       interval - Vector, 4 x 1, given as [xa xb ya yb] if the integration
%                  region is [xa, xb] x [ya, yb]
%       rule - Vector or cell, rules of Gauss quadrature of x/y, can be 
%             given as an 2 x 1 integer vector (orders of quadrature w.r.t. 
%             x and y) or a 2 x 1 cell (its values are 2-column matrices, 
%             whose columns are abscissas and weights of x/y resp.)
% 
%   [Ouput Argument]
%       output - Double, the integral
% 
% Details:
%   1. Integration region will be converted to [-1, 1]^2 by change of
%      variables before integration

function output = gq2(integrand, interval, rule)
    % Deal with sym input
    if isa(integrand, 'sym')
        integrand = matlabFunction(integrand, 'Vars', {'x','y'});
    end
    

    % Gauss quadrature rule
    if ~iscell(rule)
        [abscissa_p, weight_p] = le_gqr(rule(1));
        [abscissa_q, weight_q] = le_gqr(rule(2));
        rule = {[abscissa_p weight_p], [abscissa_q weight_q]};
    else
        [abscissa_p, weight_p] = deal(rule{1}(:, 1), rule{1}(:, 2));
        [abscissa_q, weight_q] = deal(rule{2}(:, 1), rule{2}(:, 2)); 
    end


    % Change of variables
    if ~all(interval == [-1 1 -1 1])
        range_p = interval(2)-interval(1);
        range_q = interval(4)-interval(3);
        integrand = @(p, q)(integrand((p+1)*range_p/2+interval(1), ...
                                      (q+1)*range_q/2+interval(3))...
                            *range_p*range_q/4);
    end


    % Gauss quadrature
    output = 0;
    for k = 1:size(rule{1}, 1)
        for l = 1:size(rule{2}, 1)
            output = output+integrand(abscissa_p(k), abscissa_q(l))...
                            *weight_p(k)*weight_q(l);
        end
    end

    
end