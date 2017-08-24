function [] = StratHist(hfig, chi, ID)

    if ID(6) == 'm', legstr = 'mooring'; end
    if ID(6) == 'i', legstr = 'internal'; end

    set(hfig, 'DefaultLegendBox', 'off')

    % this madness is required because h1.EdgeColor == 'auto'
    co = get(gca, 'ColorOrder');
    ico = get(gca, 'ColorOrderIndex');

    subplot(221); cla('reset')
    h1 = histogram(chi.dTdz,  'NumBins', 200, 'DisplayStyle', 'stairs', ...
                   'Normalization', 'probability', 'DisplayName', legstr, ...
                   'LineWidth', 1.5);
    hold on;
    cdf = histcounts(chi.dTdz, 'BinEdges', h1.BinEdges, ...
                     'Normalization', 'cdf');
    cbins = (h1.BinEdges(1:end-1) + h1.BinEdges(2:end))/2;
    plot(cbins, cdf, ...
         'LineStyle', '--', 'LineWidth', 1.5, 'Color', co(ico, :), ...
         'HandleVisibility', 'off')
    xlabel('dTdz')
    xlim(cbins([find_approx(cdf, 0.05) find_approx(cdf, 0.85)]))
    legend('-dynamiclegend')
    
    subplot(222)
    histogram(chi.N2,  'DisplayStyle', 'stairs', ...
              'Normalization','pdf', 'DisplayName', legstr, ...
              'LineWidth', 1.5)
    hold on
    xlabel('N2')

    subplot(223)
    histogram(log10(abs(1./chi.dTdz)),  'DisplayStyle', 'stairs', ...
              'Normalization', 'pdf', 'DisplayName', legstr, ...
              'LineWidth', 1.5)
    hold on
    xlabel('log_{10} |1/(dTdz)|')

    subplot(224)
    histogram(log10(1./chi.dTdz.^2),  'DisplayStyle', 'stairs', ...
              'Normalization', 'pdf', 'DisplayName', legstr, ...
              'LineWidth', 1.5)
    hold on
    xlabel('log_{10}1/(dTdz)^2')
end
