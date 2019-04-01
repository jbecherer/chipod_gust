function Histograms(chi, hfig, normstr, ID, legendtext)

    if ~exist('normstr', 'var'), normstr = 'pdf'; end
    if ~exist('ID', 'var'), ID = ''; end
    if ~exist('legendtext', 'var'), legendtext = ''; end
    if ~exist('hfig', 'var') | isempty(hfig), hfig = figure; end

    ID = strrep(ID, '_', '-');
    legendtext = strrep(legendtext, '_', '-');

    nbins = 300;
    chibins = linspace(-12, -1, nbins);
    epsbins = linspace(-12, 1, nbins);
    ktbins  = linspace(-8, 3, nbins);
    jqbins  = linspace(-3, 6, nbins/2);

    set(groot, 'currentfigure', hfig);
    hfig.Position(3) = 2000;
    hfig.Position(4) = 1065;
    set(hfig, 'DefaultLegendBox', 'off')

    hax(1) = subplot(221);
    set(gca, 'color', 'none')
    myhist(chi.chi, chibins, normstr, ID, legendtext)
    hold on;
    xlabel('log_{10} \chi')
    xlim([-12, -1])
    ylabel(normstr)

    hax(2) = subplot(222);
    set(gca, 'color', 'none')
    myhist(chi.eps, epsbins, normstr, ID, legendtext)
    xlim([-12, 1])
    xlabel('log_{10} \epsilon')
    ylabel(normstr)
    hold on;

    if (chi.time(11) - chi.time(10))*86400 < 2
        nt = 600;
        % 1 sec data. average Kt, Jq
        chival = moving_average(chi.chi, nt, nt);
        dTdz = moving_average(chi.dTdz, nt, nt);
        Kt = 0.5 * chival ./ dTdz;
        Jq = -1025 .* 4200 .* Kt .* dTdz;
        prefix = '10 min avg';
    else
        Kt = chi.Kt;
        Jq = chi.Jq;
        prefix = '';
    end

    hax(3) = subplot(223);
    set(gca, 'color', 'none')
    myhist(chi.Kt, ktbins, normstr, ID, legendtext)
    xlim([-8, 3])
    hold on;
    xlabel([prefix 'log_{10} K_T'])
    ylabel(normstr)

    hax(4) = subplot(224);
    set(gca, 'color', 'none')
    myhist(abs(chi.Jq), jqbins, normstr, ID, legendtext)
    xlim([-5, 6])
    hold on;
    xlabel([prefix 'log_{10} |J_q|'])
    ylabel(normstr)
end

function myhist(var, bins, normstr, ID, legendtext)

    avg = nanmean(var);
    med = nanmedian(var);

    ms = 8; % markersize

    str = sprintf('%s\n%s', legendtext, ...
                  ['\mu=' num2str(avg, '%.1e') '|mdn=' num2str(med, '%.1e')]);

    if strcmp(legendtext,'raw')   
        color = [0 0 0];
        lw = 2;
    elseif ~isempty(strfind(legendtext,'|Tz|'))
        color = [1 0 0];
        lw = 2;
    elseif ~isempty(strfind(legendtext,'n_freq'))
        color = [0 0 1];
        lw = 2;
    elseif ~isempty(strfind(legendtext, 'back flow'))
        color = [0.5 0.5 0.5];
        lw = 1;
    else
        color = choose_color(ID,'color');
        lw = choose_color(ID,'width');
    end

    if ~isempty(strfind(legendtext, 'W&DA'))
        lw = 2;
        ls = '--';
    else
        ls = '-';
    end

    if strcmpi(normstr, 'count')
        hh = histogram(log10(var), 'BinEdges', bins, 'normalization', normstr, ...
                       'displayname', str, 'displaystyle', 'stairs', ...
                       'LineWidth', lw, 'linestyle', ls, 'EdgeColor', color);
    end

    if strcmpi(normstr, 'pdf')
        hh = histogram(log10(var), 'normalization', normstr, ...
                       'displayname', str, 'displaystyle', 'stairs', ...
                       'LineWidth', lw, 'linestyle', ls, 'EdgeColor', color);
        ylim([0 max([0.5, ylim])])
    end
    hold on;

    ylims = ylim;
    plot(log10(avg), ylims(2)*0.9, 'v', 'color', hh.EdgeColor, ...
         'markersize', ms, 'handlevisibility', 'off');

    plot(log10(med), ylims(2)*0.9 , '+', 'color', hh.EdgeColor, ...
         'markersize', ms, 'handlevisibility', 'off');

    legend('-dynamiclegend', 'Location', 'northeastoutside');
end
