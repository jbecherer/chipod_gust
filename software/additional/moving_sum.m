function A = moving_sum(v,ww,ws,min_valid)
%% A = moving_sum(M,ww,ws,min_valid)
% the function calculates a moving sum (A) for the vektor (v)
% with the window width (ww)
% with the window step (ws)
% min_valid:  % of valid values in window for good average
%  length(A) = length(v)-ww

if ww == 1 | ww == 0
    A = v;
    return;
end
N = length(v);

if(N < ww)
    warning('window width is larger than vector length')
    A = nansum(v);
    return;
end

A = nan(1,round((N)/ws)-1);

if ~exist('min_valid', 'var'), min_valid = 0; end
if min_valid > 1, min_valid = min_valid/100; end

try
    AA = movsum(v, ww, 'omitnan', 'endpoints', 'discard');
    nn = movsum(~isnan(v), ww, 'omitnan', 'endpoints', 'discard');

    % decimate
    AA = AA(1:ws:end);
    nn = nn(1:ws:end);

    % at least min_valid valid values in averaging interval
    AA(nn/ww < min_valid) = NaN;
    A(1:length(AA)) = AA;
catch ME
    for i = 1:length(A)
        if((ww+(i-1)*ws)<=length(v))
            A(i) = nansum(v((1:ww)+(i-1)*ws));
        end
    end
end
