function [T1, T2] = ExtractTSFromTaoTritonPirataRama(ChipodLon, ChipodLat, ChipodDepth, ...
                                                     deployStart, deployEnd, ...
                                                     datadir,  ...
                                                     arrayName, Tres, Sres)
% This function extracts data from TAO/TRITON/RAMA/PIRATA moorings.
% Given depth of chipod it'll extract T,S from nearest two depths above &
% below the chipod.
%  [T1, T2] = ExtractFromTaoTritonPirataRama(lon, lat, depth, datadir, outputdir)
% Inputs
%            Chipod(lon,lat,depth)   -- location of chipod
%            datadir                 -- location of TaoTritonPirataRama folder
%            arrayname               -- TAO or PIRATA or RAMA
%            (Tres, Sres)            -- temporal resolution for temp / salinity.
%                                       can be '2m', '10m', '30m', 'hr', 'dy'

    ChipodDepth = abs(ChipodDepth);

    Tname = GetPmelFileName('t', ChipodLon, ChipodLat, datadir, ...
                            arrayName, Tres);
    Sname = GetPmelFileName('s', ChipodLon, ChipodLat, datadir, ...
                            arrayName, Sres);

    Tdepth = ncread(Tname, 'depth');
    Ttime = ProcessPmelTime(Tname);
    Tfull = squeeze(ncread(Tname, 'T_20'));
    Tfull(Tfull > 40) = NaN;

    Sdepth = ncread(Sname, 'depth');
    Stime = ProcessPmelTime(Sname);
    Sfull = squeeze(ncread(Sname, 'S_41'));
    Sfull(Sfull > 40) = NaN;

    % do we only have daily data?
    if strfind(Sname, 'xyzt')
        lon = ncread(Sname, 'lon');
        lat = ncread(Sname, 'lat');

        ilon = find(lon == ChipodLon);
        ilat = find(lat == ChipodLat);

        Sfull = Sfull(:,:,ilat,ilon)';
    end

    [indexT1, indexT2] = ChooseDepthLevels(Tfull(:,:), Tdepth, ChipodDepth);
    [indexS1, indexS2] = ChooseDepthLevels(Sfull(:,:), Sdepth, ChipodDepth);

    if Sdepth(indexS1) ~= Tdepth(indexT1)
        warning('Upper salinity not at same depth as temperature!');
    end
    if Sdepth(indexS2) ~= Tdepth(indexT2)
        warning('Lower salinity not at same depth as temperature!');
    end

    if isempty(deployStart), deployStart = Ttime(1); end
    if isempty(deployEnd), deployEnd = Ttime(end); end

    if Ttime(1) <= deployStart
        tind(1) = find(Ttime <= deployStart, 1, 'last');
    else
        warning('T/S data start after chipod deployment');
        tind(1) = 1;
    end

    if Ttime(end) >= deployEnd
        tind(2) = find(Ttime >= deployEnd, 1, 'first');
    else
        warning('T/S data end before chipod deployment');
        tind(2) = length(Ttime);
    end

    disp(['Chipod at ' num2str(ChipodDepth) ' m.']);
    disp(['Choosing CTDs at ' num2str(Tdepth([indexT1 indexT2])', '%dm ')]);

    T1.time = Ttime(tind(1):tind(2))';
    T1.T = Tfull(indexT1, tind(1):tind(2));
    T1.z = Tdepth(indexT1);
    T1.S = InterpAndFillSalinity(Stime, Sfull(indexS1, :), T1.time, T1.T);

    T2.time = T1.time;
    T2.T = Tfull(indexT2, tind(1):tind(2));
    T2.z = Tdepth(indexT2);
    T2.S = InterpAndFillSalinity(Stime, Sfull(indexS2, :), T2.time, T2.T);
end

function [S] = InterpAndFillSalinity(Stime, Svec, Ttime, Tvec)
    S = interp1(Stime, Svec, Ttime);

    assert(isequal(size(S), size(Tvec)));

    num_nan = sum(isnan(S) & ~isnan(Tvec));
    if num_nan < 0
        warning(['Filling ' num2str(num_nan) ' elements in Salinity']);
        S(isnan(S) & ~isnan(Tvec)) = nanmean(S);
    end
end

% find instruments with least amount of missing data within 20m of chipod
% missing data checked over chipod deployment time
function [index1, index2] = ChooseDepthLevels(var, depth, ChipodDepth)

    index1 = find( (abs(-depth + ChipodDepth) <= 10) ...
                    & (ChipodDepth > depth));
    index1 = ChooseLeastNans(var, index1);

    index2 = find( (abs(-depth + ChipodDepth) <= 10) ...
                    & (ChipodDepth < depth));
    index2 = ChooseLeastNans(var, index2);

    assert(depth(index1) < ChipodDepth & depth(index2) > ChipodDepth, ...
           'Measurements don''t span Chipod depth!');
end

function [bestIndex] = ChooseLeastNans(var, index)

    num_nans = sum(isnan(var(index,:)), 2);

    [~, bestIndex] = min(num_nans);
    bestIndex = index(bestIndex);

end