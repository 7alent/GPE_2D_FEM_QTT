% 
% Sort Character Cell
% 
% CC = CHAR_SORT(CC, REF)
%   Given a character cell cc, sort it according to reference cell ref,
%   e.g. if cc = {'d', 'b', 'f'}, ref = {'a', 'b', 'c', 'd', 'e', 'f'}, 
%   then the output will {'b', 'd', 'f'}
% 
%   [Input Argument]
%       cc - Cell, character cell to be sorted
%       ref - Cell, reference character cell
% 
%   [Ouput Argument]
%       cc - Cell, sorted character cell
% 
% Details:
%   1. This function is used in gpeq_2d()
%   2. Please ensure that cc is a subset of ref to get the correct output


function cc = char_sort(cc, ref)
    % Input check
    if ~isvector(cc) || ~isvector(ref)
        error('Inputs should be cells!');
    elseif ~iscell(cc(:)) || ~iscell(ref(:))
        error('Inputs should be cells!');
    else
        for i = cc
            if ~ischar(i{:})
                error('Inputs should be characters!');
            end
        end
        for i = ref
            if ~ischar(i{:})
                error('Inputs should be characters!');
            end
        end
    end


    % Generate indices for reference
    ref_ind = table('Size', [1 length(ref)], ...
                    'VariableTypes', repmat({'double'}, ...
                                            [1 length(ref)]), ...
                    'VariableNames', ref);
    ref_ind(1, :) = num2cell(1:length(ref));

    
    % Get indices of character cell
    cc_ind = 1:length(cc);
    for i = 1:length(cc)
        cc_ind(i) = ref_ind{1, cc{i}};
    end


    % Sort character cell
    [~, cc_ind] = sort(cc_ind);
    cc = cc(cc_ind);


end