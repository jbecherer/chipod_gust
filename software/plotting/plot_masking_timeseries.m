function plot_masking_timeseries(estimate, label, is_visible)

    if ~exist('is_visible', 'var'), is_visible = 'on'; end
    if ~exist('label', 'var'), label = ''; end

    hfig = CreateFigure(is_visible, ...
                        ['When was each mask applied? - ' strrep(label, '_', '-')]);

    [ax, axc] = create_axes(hfig, 4, 1, 0);

    ff = fieldnames(estimate.masks);

    for fff=1:length(ff)
        f = ff{fff};
        if strcmpi(f, 'dTdz') | strcmpi(f, 'N2'), axplot = ax(1); end
        if strcmpi(f, 'inst_speed') | strcmpi(f, 'back_flow'), axplot = ax(2); end
        if strcmpi(f, 'noise_floor') | strcmpi(f, 'min_n_freq'), axplot = ax(3); end
        if strcmpi(f, 'reb'), axplot = ax(3); end
        if strcmpi(f, 'ic_fit') | ~isempty(strfind(f, 'max_')), axplot = ax(4); end

        plot(axplot, estimate.time, estimate.masks.(f), ...
             'DisplayName', strrep(f, '_', '-'));
    end

    for aa = 1:length(ax)
        legend(ax(aa), '-dynamiclegend');
        axis(ax(aa), 'tight')
    end
    datetick(ax(aa), 'x', 'keeplimits')
    linkaxes(ax, 'xy')

    ylabel(ax(2), 'Number of times mask was applied for one averaged estimate')
end