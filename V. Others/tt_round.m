% 
% Rounding for TT of Multi-dimensional Arrays
%
% TT = TT_ROUND(TT, VARARGIN)
%   Conduct TT rounding process on TT representation of arbitary multi-
%   dimensional arrays tt to optimize its TT ranks
% 
%   [Input Argument]
%       tt - TT-tensor, TT-matrix or cell array, TT of the array
%       tol - Scalar, optional, tolerance of rounding (default: 1e-14)
% 
%   [Ouput Argument]
%       tt - TT-tensor, TT-matrix or cell array, rounded TT, if the input
%            tensor is TT of 1-D (vector) or 2-D (matrix), then the output 
%            will be TT-tensor/TT-matrix resp., otherwise it will be a cell
%            array
% 
% Details:
%   1. If tt is a TT-tensor or TT-matrix, this function is the same as
%      round(); while if it's not, it will be reshaped as a TT-tensor by
%      merging all mode sizes of each core and conduct TT rounding of
%      TT-tensor, then seperate merged modes
%   2. This function supposes that all mode sizes of each core are not 1,
%      otherwise this function won't get correct output
%   3. This function tends to convert tt to TT-tensor or TT-matrix and uses
%      round() first, if it fails it will get sizes of tt in cell array 
%      format


function tt = tt_round(tt, varargin)
    % Input number check
    if nargin ~= 1 && nargin ~= 2
        error('Input number should be 1 or 2!');
    end


    % Input assignment
    if nargin == 2
        tol = varargin{1};
    else
        tol = 1e-14;
    end


    % Input type check
    if iscell(tt)
        for d = 1:length(tt)
            if ~is_array(tt{d})
                error('If TT is a cell, it should be a cell array!');
            end
        end
    elseif ~isa(tt, 'tt_tensor') && ~isa(tt, 'tt_matrix')
        error('TT should be a TT-tensor, TT-matrix or cell array!');
    elseif ~is_array(tol) || ~isscalar(tol)
        error('Tolerance should be a numeric scalar!');
    elseif tol <= 0
        error('Tolerance should be positive!');
    end


    % Try to convert tensor to TT-tensor or TT-matrix
    try
        tt = tt_tensor(tt);
    catch
        try
            tt = tt_matrix(tt);
        catch
        end
    end

    
    % Rounding
    if isa(tt, 'tt_tensor') || isa(tt, 'tt_matrix')
        tt = round(tt, tol);
    else
        % Sizes of TT
        sz = tt_sz(tt);

        % Merge modes
        for d = 1:length(tt)
            tt{d} = reshape(tt{d}, ...
                            [prod(sz(d, 1:(end-2))) sz(d, (end-1):end)]);
        end

        % Squeeze first core or tt can't be converted to TT-tensor
        tt{1} = squeeze(tt{1});

        % TT-tensor rounding
        tt = tt_tensor(tt);
        tt = round(tt, tol);
        
        % Update sizes
        sz = [sz(:, 1:(end-2)) tt_sz(tt, 'rank')];

        % Seperate modes
        tt = core(tt);
        for d = 1:length(tt)
            tt{d} = reshape(tt{d}, sz(d, :));
        end
        tt{1} = squeeze(tt{1});

    end
end