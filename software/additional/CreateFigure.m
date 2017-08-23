function [hfig] = CreateFigure()

    hfig = figure('Color',[1 1 1],'visible','on', ...
                  'Position', [100 100 1400 900]);

    set(hfig, 'DefaultAxesFontSize', 18);
    set(hfig, 'renderer', 'opengl');

    set(hfig,'DefaultAxesBox','on')
    set(hfig,'DefaultAxesTickDir','out')

    set(hfig, 'DefaultHistogramEdgeColor', 'none');

end