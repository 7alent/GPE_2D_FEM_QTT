% 
% Error Plot of Finite Element Method
% 
% fem_err_plot(X, Y, BASE_OF_LOG, LABEL_X, LABEL_Y, TICKS_X, TICKS_Y, 
%              REF_LINE, TITLE_TEXT, GROUP_NAME, GROUP_INDICES, 
%              ANNOTATION_LOCATION, ANNOTATION_TEXT)
%   Plot the error of eigen values/functions in finite element method
% 
%   [Input Argument]
%       x - Vector, the data for x-axis
%       Y - Matrix, each column contains data for y-axis
%       base_of_log - Vector or scalar, base of logarithm of x and Y, to 
%                     avoid logarithm, set the 1-st/2-nd entry of it as a
%                     non-positive double or 1; if it's a scalar, base of
%                     logarithm of x and Y will be the same
%       label_x, label_Y - Character, label of x/y-axis
%       ticks_x, ticks_Y - Vector, ticks of x/y-axis
%       ref_line - Vector or matrix, for each row, the 1-st entry is the 
%                  slope (Inf for vertical line and 0 for horizontal line) 
%                  and the 2-nd one is a double in [0, 1] which controls 
%                  the intercept (if it is closer to the 1, the intercept 
%                  will be larger but the reference line will still be in 
%                  the drawing area), set it as a matrix if multiple
%                  reference lines should be drawn
%       title_text - Charactoer, text of the title
%       group_name - Character, name of the group of columns of Y
%       group_indices - Vector or cell, indices of columns of Y; indices
%                       can be double or character, set group_indices as a 
%                       character cell or double vector resp.
%       annotation_location - Vector or matrix, eahc row denotes the x-y 
%                             location of an annotation
%       annotation_text - Character or cell, text of each annotation
% 
% Details:
%   1. If base_of_log is [], x and Y data won't be processed with logarithm
%   2. If label_x/label_Y is [], ticks of x/y-axis will be generated
%      automatically
%   3. If ref_line is [], the reference line won't be drawn; for a row of
%      ref_line, if its 1-st entry is Inf, then the 2-nd entry controls the
%      abscissa of the vertical reference line (the more it is close to 1, 
%      the larger abscissa will be)
%   4. If the reference line has positive(nagetive) slope and its intercept
%      parameter (in the 2-nd column of ref_line) is set as 0/1, then the
%      bottom-right/top-left (bottom-left/top-right) points of the drawing 
%      area will be on it
%   4. If title_text is [], the title will be 'label_Y vs label_x in 
%      Different group_name'
%   5. If annotation_location/annotation_text is [], no annotation will be
%      added


function fem_err_plot(x, Y, base_of_log, label_x, label_Y, ticks_x, ...
                      ticks_Y, ref_line, title_text, group_name, ...
                      group_indices, annotation_location, annotation_text)
    % Logarithmization
    if ~isempty(base_of_log)
        if isscalar(base_of_log)
            if base_of_log > 0 && base_of_log ~= 1
                base_of_log = [base_of_log base_of_log];
            else
                base_of_log = [];
            end
        else
            if base_of_log(1) <= 0 || base_of_log(1) == 1
                base_of_log(1) = 0;
            end
            if base_of_log(2) <= 0 || base_of_log(2) == 1
                base_of_log(2) = 0;
            end
        end
    end
    if ~isempty(base_of_log)
        if base_of_log(1) > 0 && base_of_log(1) ~= 1
            x = log(x)/log(base_of_log(1));
        end
        if base_of_log(2) > 0 && base_of_log(2) ~= 1
            Y = log(Y)/log(base_of_log(2));
        end
    end


    % If Y is a vector, ensure it's a column vector
    if isrow(Y)
        Y = Y';
    end


    % Basic plotting
    marker_cell = {'pentagram', 'o', '*', 'square', 'diamond', '^', ...
                   'hexagram', 'v', '+', '>', 'x', '<', '_', '|', '.'};
    for j = 1:size(Y, 2)
        plot(x, Y(:, j), 'LineWidth', 1.5, 'Marker', ...
             marker_cell{mod(j, length(marker_cell))}, 'MarkerSize', 12)
        hold on
    end


    % Set axes labels
    set(gca, 'FontSize', 18);
    if ~isempty(base_of_log)
        if base_of_log(1) > 0 && base_of_log(1) ~= 1
            new_label_x = ['$log_{', num2str(base_of_log(1)), '}(', ...
                label_x, ')$'];
        else
            new_label_x = label_x;
        end
        if base_of_log(2) > 0 && base_of_log(2) ~= 1
            new_label_Y = ['$log_{', num2str(base_of_log(2)), '}(', ...
                label_Y, ')$'];
        else
            new_label_Y = label_Y;
        end
    else
        new_label_x = label_x;
        new_label_Y = label_Y;
    end
    xlabel(new_label_x, 'FontSize', 20, 'Interpreter', 'latex')
    ylabel(new_label_Y, 'FontSize', 20, 'Interpreter', 'latex')


    % Set axes ticks
    if isempty(ticks_x)
        x_order = 10^floor(log10(min(min(abs(x))))); % Order of magnitude
                                                     % of x
        x_digit = log10(x_order);
        x_min = (round(min(x)/10^x_digit)-1)*10^x_digit;
        x_max = (round(max(x)/10^x_digit)+1)*10^x_digit;
        xticks(x_min:(10^x_digit):x_max)
    else
        x_min = ticks_x(1);
        x_max = ticks_x(end);
        xticks(ticks_x)
    end
    if isempty(ticks_Y)
        Y_order = 10^floor(log10(min(min(abs(Y))))); % Order of magnitude 
                                                     % of Y
        Y_digit = log10(Y_order);
        Y_min = (round(min(Y)/10^Y_digit)-1)*10^Y_digit;
        Y_max = (round(max(Y)/10^Y_digit)+1)*10^Y_digit;
        yticks(Y_min:(10^Y_digit):Y_max)
    else
        Y_min = ticks_Y(1);
        Y_max = ticks_Y(end);
        yticks(ticks_Y)
    end
    axis([x_min, x_max, Y_min, Y_max])


    % Add reference line
    if ~isempty(ref_line)
        [slope, intercept] = deal(1:size(ref_line, 1));
        for i = 1:size(ref_line, 1)
            if isinf(ref_line(i, 1)) % Vertical line
                slope(i) = Inf;
                intercept(i) = (ref_line(i, 2)*(x_max-x_min)+x_min);
                ref_line_x = intercept(i)*ones(length(x), 1);
                ref_line_y = linspace(Y_min, Y_max, length(x));
            elseif ref_line(i, 1) == 0 % Horizontal line
                slope(i) = 0;
                ref_line_x = linspace(x_min, x_max, length(x));
                intercept(i) = (ref_line(i, 2)*(Y_max-Y_min)+Y_min);
                ref_line_y = intercept(i)*ones(length(x), 1);
            elseif ref_line(i, 1) > 0 % Positive slope
                slope(i) = ref_line(i, 1);
                intercept_min = Y_min-slope(i)*x_max;
                intercept_max = Y_max-slope(i)*x_min;
                intercept(i) = ref_line(i, 2)*...
                               (intercept_max-intercept_min)+intercept_min;
                ref_line_x_min = max((Y_min-intercept(i))/slope(i), x_min);
                ref_line_x_max = min((Y_max-intercept(i))/slope(i), x_max);
                ref_line_y_min = max(slope(i)*x_min+intercept(i), Y_min);
                ref_line_y_max = min(slope(i)*x_max+intercept(i), Y_max);
                ref_line_x = linspace(ref_line_x_min, ref_line_x_max, ...
                                      length(x));
                ref_line_y = linspace(ref_line_y_min, ref_line_y_max, ...
                                      length(x));
            else % Nagative slope
                slope(i) = ref_line(i, 1);
                intercept_min = Y_min-slope(i)*x_min;
                intercept_max = Y_max-slope(i)*x_max;
                intercept(i) = ref_line(i, 2)*...
                               (intercept_max-intercept_min)+intercept_min;
                ref_line_x_min = max((Y_max-intercept)/slope(i), x_min);
                ref_line_x_max = min((Y_min-intercept)/slope(i), x_max);
                ref_line_y_min = max(ref_line(i, 1)*x_max+intercept(i), ...
                                     Y_min);
                ref_line_y_max = min(ref_line(i, 1)*x_min+intercept(i), ...
                                     Y_max);
                ref_line_x = linspace(ref_line_x_min, ref_line_x_max, ...
                                      length(x));
                ref_line_y = linspace(ref_line_y_min, ref_line_y_max, ...
                                      length(x));
            end
            plot(ref_line_x, ref_line_y, 'LineWidth', 1.5, ...
                 'LineStyle', '--')
        end
    end


    % Add title
    if isempty(title_text)
        if isempty(group_name)
            title_text = ['$', label_Y, '$ vs $', label_x];
        else
            title_text = ['$', label_Y, '$ vs $', label_x, ...
                          '$ in Different ', group_name];
        end
    end
    title(title_text, 'FontSize', 24, 'Interpreter', 'latex')


    % Add legend
    if ~isempty(group_indices) && ischar(group_indices)
        group_indices = {group_indices};
    end
    if ~isempty(ref_line)
        legend_cell = cell(length(group_indices)+size(ref_line, 1), 1);
    else
        legend_cell = cell(length(group_indices), 1);
    end
    if iscell(group_indices)
        for i = 1:length(group_indices)
            legend_cell{i} = group_indices{i};
        end
    else
        for i = 1:length(group_indices)
            legend_cell{i} = {['$', group_name, ' = ', ...
                               num2str(group_indices(i)), '$']};
        end
    end
    if ~isempty(ref_line)
        for i = 1:size(ref_line, 1)
            if isinf(ref_line(i, 1))
                legend_cell{length(group_indices)+i} = ...
                    [new_label_x, ' $= ', num2str(intercept(i)), '$'];
            elseif ref_line(i, 1) == 0
                legend_cell{length(group_indices)+i} = ...
                    [new_label_Y, ' $= ', num2str(intercept(i)), '$'];
            else
                legend_cell{length(group_indices)+i} = ...
                    ['$Slope = ', num2str(ref_line(1, i)), '$'];
            end
        end
    end
    legend(legend_cell, 'FontSize', 24, 'Interpreter', 'latex', ...
           'Location', 'northeastoutside')


    % Add annotations
    if ~isempty(annotation_location) && ~isempty(annotation_text)
        if ischar(annotation_text)
            annotation_text = {annotation_text};
        end
        if size(annotation_location, 1) > 1 && ...
           size(annotation_location, 2) == 1
            annotation_location = annotation_location';
        end
        for i = 1:length(annotation_text)
            text(annotation_location(i, 1), annotation_location(i, 2), ...
                 annotation_text{i}, 'FontSize', 20, 'Interpreter', ...
                 'latex')
        end
    end

    
    hold off


end