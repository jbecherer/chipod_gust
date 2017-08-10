function [v0] = fit_pitot_v0(a_time, a_spd, p_time, p_v_cal, s_pd, varargin)
%% [v0] = fit_pitot_v0(a_time, a_spd, p_time, p_v_cal, s_pd, [do_plot])
%     
%     This function calculates v0 by fiting calibrated pitot voltages against
%     a given velocity time series
%
%     INPUT
%        a_time      :  time vector of reference speed
%        a_spd       :  reference speed vector 
%        p_time      :  time vector of pitot 
%        p_v_cal     :  pitot voltage (note this voltage should be calibrated
%                       for P, T, and tilt already!)
%        s_pd        :  calibration slope from voltage into dynamic pressure [Pa/volt]
%        do_plot     :  optional argument to produce a figure (default do_plot=0)
%
%     OUTPUT
%        v0          :  matching v0 for P_d = s_pd*(v_cal-v0)
%
%   created by: 
%        Johannes Becherer
%        Tue Nov  1 10:32:38 PDT 2016

if nargin == 6
   do_plot = varargin{1};
else
   do_plot = 0;
end
%---------------------intrepolate referece speed to pitot time----------------------

p_ref_spd = interp1(a_time, a_spd, p_time);

%---------------------cal v0 time series----------------------
tmp_v0 = p_v_cal - .5*1025/s_pd*p_ref_spd.^2;

%---------------------average for a single value----------------------
v0 = nanmean(tmp_v0);

%---------------------plot if requested----------------------
if do_plot
    CreateFigure;
    set(gcf, 'DefaultLineLineWidth', 1);
    ax(1) = subplot(211);
      plot(p_time, tmp_v0);
         hold all;
         xl = get(gca,'Xlim');
      plot(xl, [1 1]*v0);
      legend('pitot volt - ref. spd. volt', 'inferred V_0')
      xlim([nanmin(p_time) nanmax(p_time)])
      datetick('keeplimits')
      title('determine V_0 by fitting to ADCP')

   ax(2) = subplot(212);
   hold on;
   plot(p_time, .5*1025/s_pd*p_ref_spd.^2);
   plot(p_time, p_v_cal - v0);
   plot(get(gca, 'XLim'), [0 0], '--', 'color', [1 1 1]*0.6)
   xlim([nanmin(p_time) nanmax(p_time)])
   title(['num(negative speeds) = ' num2str(sum(p_v_cal-v0 < 0))])
   legend('pitot calibrated voltage', 'ref spd converted to voltage', ...
          'pitot calibrated voltage with offset correction')
   datetick('keeplimits')

   linkaxes(ax, 'x')
end

