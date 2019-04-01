function [] = StratHist(hfig, chi, ID)

    if ID(2) == 'm', legstr = 'mooring'; end
    if ID(2) == 'i', legstr = 'internal'; end
    if ID(2) == 'w', legstr = 'sorted'; end

    set(groot, 'currentfigure', hfig);
    set(hfig, 'DefaultLegendBox', 'off')

    % this madness is required because h1.EdgeColor == 'auto'
    co = get(gca, 'ColorOrder');
    ico = get(gca, 'ColorOrderIndex');

    subplot(121);
    h1 = histogram(chi.dTdz,  'NumBins', 200, 'DisplayStyle', 'stairs', ...
                   'Normalization', 'probability', 'DisplayName', legstr, ...
                   'LineWidth', 1.5);
    hold on;
    % cdf = histcounts(chi.dTdz, 'BinEdges', h1.BinEdges, ...
    %                  'Normalization', 'cdf');
    % cbins = (h1.BinEdges(1:end-1) + h1.BinEdges(2:end))/2;
    % plot(cbins, cdf, ...
    %      'LineStyle', '--', 'LineWidth', 1.5, 'Color', co(ico, :), ...
    %      'HandleVisibility', 'off')
    xlabel('dTdz')
    ylabel('pdf');
    % xlim(cbins([find_approx(cdf, 0.05) find_approx(cdf, 0.85)]))
    legend('-dynamiclegend')

    subplot(122)
    histogram(2*pi./sqrt(chi.N2(chi.N2 > 0))/60,  'BinEdges', [0:1:90], ...
              'DisplayStyle', 'stairs', ...
              'Normalization','pdf', 'DisplayName', legstr, ...
              'LineWidth', 1.5)
    hold on
    xlabel('Buoyancy period (2*\pi/N min)')
    ylabel('pdf')
    xlim([0 90])
    set(gca, 'XTick', [0:5:90])
    legend('-dynamiclegend')

    % subplot(223)
    % histogram(log10(abs(1./chi.dTdz)),  'DisplayStyle', 'stairs', ...
    %           'Normalization', 'pdf', 'DisplayName', legstr, ...
    %           'LineWidth', 1.5)
    % hold on
    % xlabel('log_{10} |1/(dTdz)|')
    % ylabel('pdf')
    % legend('-dynamiclegend')

    % subplot(224)
    % histogram(log10(1./chi.dTdz.^2),  'DisplayStyle', 'stairs', ...
    %           'Normalization', 'pdf', 'DisplayName', legstr, ...
    %           'LineWidth', 1.5)
    % hold on
    % xlabel('log_{10}1/(dTdz)^2')
    % ylabel('pdf')
end
