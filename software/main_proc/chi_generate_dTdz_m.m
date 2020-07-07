function [Tz_m] = chi_generate_dTdz_m(t1, z1, T1, S1, t2, z2, T2, ...
                                      S2,  sdir, use_TS_relation, ...
                                      ChipodDepth, vis)
% [Tz_m] = chi_generate_dTdz_m(t1, z1, T1, S1, t2, z2, T2, S2,  sdir);
%
%        This function generates an input file for chi-processing dTdz_m.m at
%        directory sdir based on 2 CTD-time series T1 above chipod/gust and T2
%        below
%        
%        Input
%           t1          : time vector of 1st CTD time series
%           z1          : depth in meter from surface
%           T1          : temperature vector of 1st CTD time series
%           S1          : salinity vector of 1st CTD time series (could also be constant scalar)
%           t2          : time vector of 2nd CTD time series
%           z2          : depth in meter from surface
%           T2          : temperature vector of 2nd CTD time series
%           S2          : salinity vector of 2nd CTD time series (could also be constant scalar)
%           sdir        : directory to save dTdz_m.met to
%       use_TS_relation : if 1, fit TS relation; else do naive differencing
%          vis          :  visible figure 'on' 'off'
% 
%        Output
%           Tz_m.time   : time vector  (1min averages)
%           Tz_m.Tz     : temperature gradient
%           Tz_m.N2     : buoyancy freqency
%
%   created by: 
%        Johannes Becherer
%        Fri Sep  2 13:49:20 PDT 2016


% make figure unvisisible
if nargin < 12
   vis = 'on';
end

picdir   =  [sdir '../pics/'];



%---------------------check salinity id sclar----------------------
   if length(S1)==1
      S1 = ones(size(T1))*S1;
   end
   if length(S2)==1
      S2 = ones(size(T2))*S2;
   end

   if ~exist('use_TS_relation', 'var')
       use_TS_relation = 0;
   end

%---------------------create time vector----------------------
   % find beginning
   ts = max(t1(1), t2(1));
   % find end
   tf = min(t1(end), t2(end));

   % construct time array
   dt = 1/(24*60);      % one minute
   time = ts:dt:tf;     % time array

   % find delta time for the CTD series
   dt1 = nanmean(diff(t1));
   dt2 = nanmean(diff(t2));

   % vertical distance of both sensors in meter
   dz = nanmean(abs(z1-z2));

%---------------------cal gradient----------------------
   % in case the time series are sampled much quicker then one minute 
   % low pass filter them
   if dt1/dt < 1
      T1 = qbutter(T1, dt1/dt);
      S1 = qbutter(S1, dt1/dt);
   end
   if dt2/dt < 1
      T2 = qbutter(T2, dt2/dt);
      S2 = qbutter(S2, dt2/dt);
   end

   % interpolate time series on common time stamp
   T1_int = interp1(t1, T1, time);
   T2_int = interp1(t2, T2, time);
   S1_int = interp1(t1, S1, time);
   S2_int = interp1(t2, S2, time);

   % cal temperature gradient
   Tz_m.Tz = (T1_int-T2_int)/dz;
   % cal salinity gradient
   Tz_m.Sz = (S1_int-S2_int)/dz;

   if ~use_TS_relation
       if length(z1) ~= 1
           z1_int = interp1(t1, z1, time);
       else
           z1_int = z1
       end
       if length(z2) ~= 1
           z2_int = interp1(t2, z2, time);
       else
           z2_int = z2
       end

       % cal density
       D1 = sw_pden(S1_int,T1_int,abs(z1_int), 0.5*abs(z1_int+z2_int));
       D2 = sw_pden(S2_int,T2_int,abs(z2_int), 0.5*abs(z1_int+z2_int));

       % N2
       g        = 9.81;
       rho_0    = 1025;
       Tz_m.N2  = -g/rho_0*(D1-D2)/dz;
   else
       [N2_1,Sz_1,~]  = cal_N2_from_TS(t1, T1,  S1, ones(size(S1))*abs(.5*(z1+z2)), time, Tz_m.Tz, 7200, ChipodDepth);
       [N2_2,Sz_2,~]  = cal_N2_from_TS(t2, T2,  S2, ones(size(S1))*abs(.5*(z1+z2)), time, Tz_m.Tz, 7200, ChipodDepth);
       Tz_m.N2 = (N2_1 + N2_2)/2;
       Tz_m.Sz = (Sz_1 + Sz_2)/2;
   end

   Tz_m.time  = time;

   numneg = length(Tz_m.N2(Tz_m.N2 < 0));
   % some plotting
    fig = figure('Color',[1 1 1],'visible',vis,'Paperunits','centimeters',...
            'Papersize',[30 20],'PaperPosition',[0 0 30 20]);
    plot(Tz_m.time, Tz_m.N2, 'b'); hold on;
    plot(Tz_m.time(Tz_m.N2 < 0), Tz_m.N2(Tz_m.N2 < 0), 'r.');
    plot(xlim, [0, 0], '--', 'color', [1 1 1]*0.6);
    ylabel('N^2')
    title(['Negative N^2 in red = ' num2str(numneg) '/' num2str(length(Tz_m.N2)) ...
           ' = ' num2str(numneg/length(Tz_m.N2) * 100) '%'])
    datetick

    if exist(picdir)
         print(fig,[picdir 'Tz_m.png'],'-dpng','-r200','-painters')
      end

%---------------------save----------------------
   save([sdir 'dTdz_m.mat'], 'Tz_m');
