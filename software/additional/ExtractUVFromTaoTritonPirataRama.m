function [moor] = ExtractUVFromTaoTritonPirataRama(ChipodLon, ChipodLat, ChipodDepth, ...
                                                   deployStart, ...
                                                   deployEnd, datadir, ...
                                                   arrayName, res)
% This function extracts data from TAO/TRITON/RAMA/PIRATA moorings.
% Given depth of chipod it'll extract U,V from nearest two depths above &
% below the chipod.
%  [T1, T2] = ExtractFromTaoTritonPirataRama(lon, lat, depth, datadir, outputdir)
% Inputs
%            Chipod(lon,lat,depth)   -- location of chipod
%            datadir                 -- location of TaoTritonPirataRama folder
%            arrayname               -- TAO or PIRATA or RAMA
%            res                     -- temporal resolution.
%                                       can be '2m', '10m', '30m', 'hr', 'dy'


    ChipodDepth = abs(ChipodDepth);

    fname = GetPmelFileName('cur', ChipodLon, ChipodLat, datadir, ...
                            arrayName, res);
    moor.depth = ncread(fname, 'depth');
    assert(length(moor.depth) == 1, ...
           ['Velocity available at more than one depth. Giving up!']);
    moor.time = ProcessPmelTime(fname);

    if isempty(deployStart), deployStart = Ttime(1); end
    if isempty(deployEnd), deployEnd = Ttime(end); end

    if moor.time(1) <= deployStart
        tind(1) = find(moor.time < deployStart, 1, 'last');
    else
        warning('Current data start after chipod deployment');
        tind(1) = 1;
    end

    if moor.time(end) >= deployEnd
        tind(2) = find(moor.time > deployEnd, 1, 'first');
    else
        warning('Current data end before chipod deployment');
        tind(2) = length(moor.time);
    end

    if tind(1) == tind(2)
        error('No velocity data during deployment!');
    end

    moor.u = squeeze(ncread(fname, 'U_320'))/100; moor.u(abs(moor.u) > 10) = NaN;
    moor.v = squeeze(ncread(fname, 'V_321'))/100; moor.v(abs(moor.v) > 10) = NaN;

    disp(['Chipod at ' num2str(ChipodDepth) ' m.']);
    disp(['Choosing velocity at ' num2str(moor.depth, '%dm ')]);

    % since I have only one depth
    moor.depth = ChipodDepth;

    moor.u = moor.u(tind(1):tind(2));
    moor.v = moor.v(tind(1):tind(2));
    moor.time = moor.time(tind(1):tind(2));
end