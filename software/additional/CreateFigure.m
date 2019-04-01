function [hfig] = CreateFigure(vis, name)
    if nargin <1
        vis = 'on';
    end

    hfig = figure('Color',[1 1 1],'visible', vis, ...
                  'Position', [100 100 1400 900]);

    set(hfig, 'DefaultAxesFontSize', 16);
    set(hfig, 'renderer', 'opengl');

    set(hfig,'DefaultAxesBox','on')
    set(hfig,'DefaultAxesTickDir','out')

    set(hfig, 'DefaultHistogramEdgeColor', 'none');

    if exist('name', 'var')
        hfig.Name = name;
        annotation('textbox', [0.1, 0.9, 0.1, 0.1], ...
                   'String', strrep(name, '_', '-'), ...
                   'FitBoxToText', 'on', 'HorizontalAlignment', 'left', ...
                   'FontWeight', 'bold', 'LineStyle', 'none', 'FontSize', 16);
    end
end