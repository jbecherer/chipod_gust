% [ax] = plot_estimate(chi, name, window, hfig, t0, t1)
%
% Quick way to plot an estimate
% Inputs:
%        chi : structure containing estimate
%        name : legend label for estimate
%        window : optional, averaging window (seconds), none by default
%        hfig: optional, figure handle, calls gcf() if not provided
%        t0, t1 : optional, time subset the plot
%
% Outputs:
%        ax : axes handles

function [ax] = plot_estimate(chi, name, window, hfig, t0, t1)

    if ~exist('name', 'var'), name = 'chi'; end
    if ~exist('window', 'var'), window = 0; end
    if ~exist('t0', 'var') | isempty(t0), t0 = nanmin(chi.time); end
    if ~exist('t1', 'var') | isempty(t1), t1 = nanmax(chi.time); end
    if ~exist('hfig', 'var') | isempty(hfig), hfig = gcf(); end

    i0 = find_approx(chi.time, t0, 1);
    i1 = find_approx(chi.time, t1, 1);
    tind = i0:i1;

    dt = (chi.time(2) - chi.time(1))*86400;
    ww = round(window/dt);
    time = moving_average(chi.time(tind), ww, ww);

    try
        set(groot, 'currentfigure', hfig);
        ax = findobj(hfig.Children, 'type', 'axes');
    catch ME
        hfig = CreateFigure;
    end

    if isempty(ax)
        % clf(hfig);
        [ax, ~] = create_axes(hfig, 5, 1, 0);
    end

    set(hfig, 'currentaxes', ax(1))
    semilogy(time, moving_average(chi.chi(tind), ww, ww), 'displayname', name)
    ylabel('\chi')
    set(ax(1), 'yscale', 'log');
    Common()
    ylim([1e-12, 1e-3]);
    grid on;

    set(hfig, 'currentaxes', ax(2))
    old_ylim = ylim;
    Tzavg = moving_average(chi.dTdz(tind), ww, ww);
    plot(time, Tzavg, 'displayname', name)
    new_ylim = [0.98*min(Tzavg), 1.05*max(Tzavg)];
    if all(isnan(new_ylim)), error('plot_estimate: dTdz All NaNs'); end
    if diff(new_ylim) < 0  % for negative limits
        new_ylim = fliplr(new_ylim);
    end
    if isequal(old_ylim, [0, 1])
        ylim(new_ylim);
    else
        ylim([min([old_ylim, new_ylim]), max([old_ylim, new_ylim])])
    end
    hold on;
    plot(xlim, [0, 0], 'k--');
    ylabel('dT/dz')
    Common()
    %symlog(gca, 'y', 5e-3);

    set(hfig, 'currentaxes', ax(3))
    try
        semilogy(time, moving_average(chi.eps(tind), ww, ww), 'displayname', name)
    catch ME
        semilogy(time, moving_average(chi.eps1(tind), ww, ww), 'displayname', name)
    end
    ylabel('\epsilon')
    set(ax(3), 'yscale', 'log');
    ylim([10.^[-12, -3]])
    grid on;
    Common()

    set(hfig, 'currentaxes', ax(4))
    if isfield(chi, 'Kt')
        semilogy(time, moving_average(chi.Kt(tind), ww, ww), ...
                 'displayname', name);
    end
    ylabel('K_t')
    set(ax(4), 'yscale', 'log');
    ylim([10.^[-7.5, 0]])
    grid on;
    Common()

    set(hfig, 'currentaxes', ax(5))
    if isfield(chi, 'Jq')
        plot(time, ...
             moving_average(chi.Jq(tind), ww, ww), 'displayname', name)
    end
    ylabel('J_q^t')
    Common()
    legend('-DynamicLegend');

    linkaxes(ax, 'x')
    if t0 < t1
        xlim([t0, t1])
    end
    datetick('x', 'mmm-dd HH:MM', 'keeplimits')

    for aa=1:length(ax)-1
        set(ax(aa), 'xticklabel', []);
    end
end

function Common()
    hold on
    datetick('x', 'mm-dd HH:MM', 'keeplimits')
end