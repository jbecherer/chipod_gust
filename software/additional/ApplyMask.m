function [chi] = ApplyMask(chi, maskvar, relation, criterion, maskname, Trange)

    if ~exist('Trange', 'var'), Trange = 1:length(chi.chi); end

    if relation == '>'
        mask = maskvar > criterion;
    elseif relation == '<'
        mask = maskvar < criterion;
    end

    numnans = sum(isnan(chi.chi(Trange)));
    chi.chi(mask) = nan;
    chi.eps(mask) = nan;
    chi.Kt(mask) = nan;
    chi.Jq(mask) = nan;
    chi.mask = isnan(chi.chi);

    numnewnans = sum(isnan(chi.chi(Trange)));

    disp([maskname ' ' relation ' ' num2str(criterion, '%1.1e') ...
          ' NaN-ed out ' num2str((numnewnans-numnans)/ ...
                                 length(Trange)*100, '%.2f') ...
          '% of estimates'])

end
