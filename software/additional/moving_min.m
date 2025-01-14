function A = moving_min(v,ww,ws)
%% A = moving_min(M,ww,ws)
% the function calculates a moving minimum (A) for the vektor (v)
% with the window width (ww)
% with the window step (ws)
%  length(A) = length(v)-ww

if ww == 1 | ww == 0
    A = v;
    return;
end
N = length(v);

if(N < ww)
    error('window width is larger than vector length')
end

A = nan(1,round((N)/ws)-1);

try
    AA = movmin(v, ww, 'omitnan', 'endpoints', 'discard');
    AA = AA(1:ws:end);
    A(1:length(AA)) = AA;
catch ME:
    for i = 1:length(A)
      if((ww+(i-1)*ws)<=length(v))
        A(i) = nanmin(v((1:ww)+(i-1)*ws)); 
      end
    end
end

end
