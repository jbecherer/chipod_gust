function [] = TestMask(chi, mask_var, mask_sign, mask_vals, name, t0, t1)

    if exist('t0', 'var')
        ff = fieldnames(chi);
        i0 = find_approx(chi.time, t0);
        i1 = find_approx(chi.time, t1);

        for f=1:length(ff)
            chi.(ff{f}) = chi.(ff{f})(i0:i1);
        end
    end

    chiold = chi;

    hfig = CreateFigure;
    Histograms(chi, hfig, 'count', 'raw');
    for ii=1:length(mask_vals)
        chi = ApplyMask(chi, mask_var, mask_sign, mask_vals(ii), name);
        if mask_sign == '>'
            label = [name ' < ' num2str(mask_vals(ii), '%.1e')];
        else
            label = [name ' > ' num2str(mask_vals(ii), '%.1e')];
        end
        Histograms(chi, hfig, 'count', label);

        chi = chiold; % reset
    end

    legend('-dynamiclegend')