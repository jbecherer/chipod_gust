% Given appropriate temperature gradients in wda, estimate volume avg Kt, Jq
%        [wda] = process_wda_estimate(chi, chi_wda)
%
% Inputs:
%        chi : nominally chi structure to process chi.
%              Alternatively can provide 1 sec averaged T to estimate just the ...
%              mean sorted gradient
%              At minimum must contain 'T' and 'time' subfields.
%
% Outputs:
%        wda : Structure with subfields 'time', 'dTdz' at least.
%              If provided with chi, will return 'Jq', 'Kt', 'eps', 'chi' subfields

function [wda] = process_wda_estimate(chi, chi_wda)

    % This is because you might provide 1 sec T as input
    % to pre-process sorted gradient
    chi_is_empty = ~isfield(chi, 'chi');

    ndt = round(chi_wda.dt(1) / diff(chi.time(1:2)*86400));

    if ~chi_is_empty
        wda.Jq = nan(1, length(chi_wda.tstart));
        wda.Kt = nan(1, length(chi_wda.tstart));
        wda.eps = nan(1, length(chi_wda.tstart));
        wda.chi = nan(1, length(chi_wda.tstart));
    end

    wda.dTdz = nan(1, length(chi_wda.tstart));
    wda.dzdT = nan(1, length(chi_wda.tstart));
    if isfield(chi_wda, 'no_min_dz')
        chi_wda.no_min_dz(isnan(chi_wda.no_min_dz)) = 0;
        wda.no_min_dz = chi_wda.no_min_dz';
    end

    % at this point, chi might have been trimmed (it is trimmed in
    % combine_turbulence); so we need to find the right time chunks in wda so
    % that things line up right.
    % chunk0 is the wda index of the time chunk that is closest to the first good chi observation
    chunk0 = find(chi_wda.tstart >= chi.time(1), 1, 'first');
    chit0 = find_approx(chi.time, chi_wda.tstart(chunk0));

    nnan = 0;
    for tt=chunk0:length(chi_wda.tstart)
        maxt = min(chit0+2*ndt, length(chi.time));

        % narrow the search space for speed using ndt
        chit0 = find_approx(chi.time(chit0:maxt), chi_wda.tstart(tt)) + chit0-1;
        chit1 = find_approx(chi.time(chit0+1:maxt), chi_wda.tstop(tt), 2) + chit0;

        % I can only start exiting early here because this loop depends on
        % chit0 from the previous iteration
        if isfield(wda, 'no_min_dz') & wda.no_min_dz(tt), continue; end

        % if at the end, they probably don't line up; discard this one point
        if chit1(1) == length(chi.time), break; end

        if abs(chi.time(chit1(1)) - chi_wda.tstop(tt))*86400 > 0.5
            chit1 = chit1(2);
        else
            chit1 = chit1(1);
        end

        % make sure the time chunks are lining up
        % disp(['Compare chit0 = ' datestr(chi.time(chit0)) ' vs tstart = ' datestr(chi_wda.tstart(tt))]);
        % disp(['Compare chit1 = ' datestr(chi.time(chit1)) ' vs tstop = ' datestr(chi_wda.tstop(tt))]);
        if (abs(chi.time(chit0) - chi_wda.tstart(tt)) > 0.5/86400) ...
           | (abs(chi.time(chit1) - chi_wda.tstop(tt)) > 0.5/86400)
            continue;
        end

        if ~chi_is_empty
            chisub = chi.chi(chit0:chit1);

            chi_is_nan = isnan(chisub);
            if all(chi_is_nan), continue; end

            % if we have chi estimates for less than 30% of the time interval,
            % note that and carry along, result is NaN
            if sum(~chi_is_nan)/length(chisub) < 0.3
                % if the valid estimates are all 0, then we keep that.
                if ~all(chisub(~chi_is_nan) == 0)
                    nnan = nnan + 1;
                    continue;
                end
            end
        end

        Tchi = chi.T(chit0:chit1);

        Tbins = chi_wda.Tbins(:, tt);
        Tbins = Tbins(~isnan(Tbins));

        if length(Tbins) > 1
            dz = chi_wda.dz(1:length(Tbins)-1, tt);
            dTdz = chi_wda.dTdz_bins(1:length(Tbins)-1, tt);
            wda.dTdz(tt) = sum(dz .* dTdz, 'omitnan')./sum(dz, 'omitnan');
            wda.dzdT(tt) = 1./(sum(dz .* 1./dTdz, 'omitnan')./sum(dz, 'omitnan'));

            if ~chi_is_empty
                avg = isoscalar_average([chi.chi(chit0:chit1); chi.eps(chit0:chit1)], ...
                                        Tchi, Tbins);
                chiavg = avg(1, :);
                epsavg = avg(2, :);

                Jqvec = -chiavg' ./ dTdz * 4200 * 1025 * 0.5;
                Ktvec = chiavg' ./ dTdz.^2 * 0.5;

                wda.chi(tt) = sum(dz .* chiavg', 'omitnan')./sum(dz, 'omitnan');
                wda.eps(tt) = sum(dz .* epsavg', 'omitnan')./sum(dz, 'omitnan');
                wda.Jq(tt) = sum(dz .* Jqvec, 'omitnan')./sum(dz, 'omitnan');
                wda.Kt(tt) = sum(dz .* Ktvec, 'omitnan')./sum(dz, 'omitnan');
            end
        end
    end

    wda.time = (chi_wda.tstart + chi_wda.tstop)/2;
    if ~chi_is_empty
        disp(['      WDA: ' num2str(nnan) ' = ' num2str(nnan/length(wda.chi)*100) ...
              '% time intervals have chi=NaN. '])
        disp(['      WDA: ' num2str(sum(wda.no_min_dz)) ' = ' ...
              num2str(sum(wda.no_min_dz)/length(wda.chi)*100) ...
              '% time intervals did not see enough pumping to make an ' ...
              'estimate. ']);

        wda.N2 = interp1(chi.time(~isnan(chi.time)), chi.N2(~isnan(chi.time)), wda.time);
        wda.T = interp1(chi.time(~isnan(chi.time)), chi.T(~isnan(chi.time)), wda.time);
        wda.eps_Kt = wda.N2 .* wda.Kt/0.2;
        wda.eps_Kt(wda.N2 < 0) = nan;
    end
end