function [] = TestMask(chi, mask_var, mask_sign, mask_vals, name)

    hfig = CreateFigure;
    chiold = chi;
    Histograms(chi, hfig, 'count', 'raw');
    for ii=1:length(mask_vals)
        chi = ApplyMask(chi, mask_var, mask_sign, mask_vals(ii), name);
        if mask_sign == '>'
            label = [name ' < ' num2str(mask_vals(ii), '%.1e')];
        elseif mask_sign == '<'
            label = [name ' > ' num2str(mask_vals(ii), '%.1e')];
        elseif mask_sign == '='
            label = [name ' ~= ' num2str(mask_vals(ii), '%.1e')];
        end
        Histograms(chi, hfig, 'count', label);

        chi = chiold; % reset
    end

    legend('-dynamiclegend')