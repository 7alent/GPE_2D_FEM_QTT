% 
% Get All Sizes of Cores of a TT
% 
% SZ = TT_SZ(TT)
%   Get the sizes of all the TT cores of a TT, the output will be a matrix
%   whose rows correspond to cores and the last 2 columns denote TT ranks
%   Supposed all mode sizes of cores are not less than 2, while rank sizes
%   can be 1
% 
%   [Input Argument]
%       tt - TT-tensor, TT- matrix or cell array, the TT
%       sz_type - Character, optional, if given, it should be either 'mode'
%                 or 'rank', the output will only contain mode/rank sizes
%                 resp., if not given, output will contain all sizes 
%                 (default: [])
% 
%   [Ouput Argument]
%       sz - Matrix, sizes of the cores
% 
% Details:
%   1. This function tends to convert tt to TT-tensor or TT-matrix and uses
%      round() first, if it fails it will get sizes of tt in cell array 
%      format


function sz = tt_sz(tt, varargin)
    % Input number check
    if nargin ~= 1 && nargin ~= 2
        error('Input number should be 1 or 2!');
    end


    % Input assignment
    if nargin == 2
        sz_type = varargin{1};
    else
        sz_type = [];
    end


    % Input type check of tensor
    if iscell(tt)
        for i = 1:length(tt)
            if ~is_array(tt{i})
                error(['If tensor input is a cell, it should be a ', ...
                       'cell array!']);
            end
        end
    elseif ~isa(tt, 'tt_tensor') && ~isa(tt, 'tt_matrix')
        error('Wrong type of tensor input!');
    end
    

    % Input type check of size type
    if ~isempty(sz_type)
        if ~ischar(sz_type)
            error('Size type should be a character!');
        elseif ~strcmp(sz_type, 'mode') && ~strcmp(sz_type, 'rank')
            error('Size type should be ''mode'' or ''rank''!');
        end
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


    % Get sizes
    if iscell(tt)
        d = length(tt); % NUmber of cores
        m = length(size(tt{1}))+1; % Number of mode and rank sizes of each 
                                   % core
        sz = zeros(d, m);
        sz(1, m-1) = 1;
        sz(d, m) = 1;
        for i = 1:length(tt)
            if i == 1
                sz(i, [1:(m-2) m]) = size(tt{i});
            elseif i == d
                sz(i, 1:(m-1)) = size(tt{i});
            else
                sz(i, :) = size(tt{i});
            end
        end
    elseif isa(tt, 'tt_tensor')
        m = 3;
        sz = [tt.n tt.r(1:tt.d) tt.r(2:(tt.d+1))];
    else
        m = 4;
        sz = [tt.n tt.m tt.r(1:tt.d) tt.r(2:(tt.d+1))];
    end


    % Select certain type of sizes
    if ~isempty(sz_type)
        if strcmp(sz_type, 'mode')
            sz = sz(:, 1:(m-2));
        else
            sz = sz(:, (m-1):m);
        end
    end

    
end