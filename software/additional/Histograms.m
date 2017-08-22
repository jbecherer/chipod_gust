function Histograms(chi, hfig, normstr, legstr)

    figure(hfig)
    hfig.Position(3) = 2000;
    set(hfig, 'DefaultLegendBox', 'off')

    subplot(221);
    set(gca, 'color', 'none')
    myhist(chi.chi, normstr, legstr)
    hold on;
    xlabel('log_{10} \chi')
    ylabel(normstr)

    subplot(222);
    set(gca, 'color', 'none')
    myhist(chi.eps, normstr, legstr)
    xlabel('log_{10} \epsilon')
    ylabel(normstr)
    hold on;

    subplot(223)
    set(gca, 'color', 'none')
    myhist(chi.Kt, normstr, legstr)
    hold on;
    xlabel('log_{10} K_T')
    ylabel(normstr)

    subplot(224)
    set(gca, 'color', 'none')
    myhist(abs(chi.Jq), normstr, legstr)
    hold on;
    xlabel('log_{10} |J_q|')
    ylabel(normstr)
end

function myhist(var, normstr, legstr)
    nbins = ceil(sqrt(numel(var)));
    legstr = getstats(var, legstr);
    histogram(log10(var), nbins, 'normalization', normstr, ...
              'displayname', legstr, ...
              'displaystyle', 'stairs', 'LineWidth', 1.5);
    legend('-dynamiclegend');
end

function str = getstats(var, strin)

    str = [' | \mu=' num2str(nanmean(var), '%.1e')];
    str = [str ', mdn=' num2str(nanmedian(var), '%.1e')];
    str = [strin str];
end