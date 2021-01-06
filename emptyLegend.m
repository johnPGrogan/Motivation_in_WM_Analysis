function h = emptyLegend(n, opts, cols, names, legOpts)
% function h = emptyLegend(n, opts, cols, names)
% plots NaN so that you have plot handles for an empty plot which can be
% used for legends
% 
% Inputs:
%   n = number of lines/plots to make
%   opts = [N x 1] cell array of strings to pass to plot e.g. '-x'
%   cols = [N x 1] cell array of colour args to pass to plot(NaN,'Color'...)
%   names = cell array of legend names
%   legOpts = cell array of options passed to legend


    if ~exist('legOpts','var') || isempty(legOpts)
        legOpts = {};
    end
    
    % plot each set of options as NaN
    for i = 1:n
        h(i) = plot(NaN, opts{i}, 'Color', cols{i}, 'LineWidth',2, 'MarkerFaceColor', cols{i});
        hold on;
    end
    
    legend(h, names(1:n), legOpts{:}); % make legend


end