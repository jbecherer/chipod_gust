%     this script is meant to combine all processed turbulence information 
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________set flags______________________
   do_combine  =  1; % do actually the averaging calculation (can take a couple o minutes)
   do_plot     =  1; % generate a comparison plot between the different estimates 
   do_mask     =  1; % NaN chi estimates using min dTdz, speed thresholds

   % set thresholds for masking
   min_N2 = 1e-9;
   min_dTdz = 1e-4;
   min_spd = 0.05;
   min_inst_spd = min_spd; % min instantaneous speed past sensor
   mask_dTdz = ''; % '' for none, 'm' for mooring, 'i' for internal
                    % in addition to whats in chi.dTdz
   mask_inst_spd = 1; % estimates are crappy if sensor isn't moving
                      % enough.
                      % screws the spectrum calculation...
   mask_spd = ''; % masking such that background flow flushes
                  % sensed water volume
                  % 'm' for mooring, 'p' for pitot,
                  % '' to choose based on what was used in chi estimate
   avgwindow = 600; % averaging window in seconds

   ChipodDepth = 30;

   normstr = 'count'; % normalization for histograms

   % if you want to restrict the time range that should be combined
   % use the following
   % This restricts time range for BOTH sensors on chipods
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

   % if one sensor dies earlier, specify that time here.
   T1death = datenum(2060, 1, 1, 0, 0, 0); % chipod T1 or gustT T sensor
   T2death = datenum(2060, 1, 1, 0, 0, 0); % T2 sensor

   % additional time ranges to NaN out as necessary
   % make an array that looks like
   % nantimes{sensor_number} = [start_time1, end_time1;
   %                            start_time2, end_time2;]
   % start_time & end_time must be datenum
   nantimes{1} = []; % sensor T1 for chipod or T sensor on gusT
   nantimes{2} = []; % sensor T2 for chipod
   nantimes{3} = []; % pitot sensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

   % determine if chipod or gusT
   load ../calib/header.mat
   if isfield(head.coef, 'T1')
       isChipod = 1;
   else
       isChipod = 0;
   end

if do_mask
    if mask_dTdz == 'm'
        disp('additional masking using mooring dTdz')
        load([basedir 'input/dTdz_m.mat'])
        Tz = Tz_m;
        clear Tz_m;

    elseif mask_dTdz == 'i'
        disp('additional masking using internal dTdz')
        load([basedir 'input/dTdz_i.mat'])
        Tz = Tz_i;
        clear Tz_i;
    else
        Tz = [];
    end

    mask_spd_initial = mask_spd;
end

%_____________________find all available chi data______________________
if(do_combine)
   dirname = [savedir 'chi/'];
   d = dir(dirname);
   % eliminate directories
   d = d(~[d(:).isdir]);

   % run through all chi files
   for i = 1:length(d)

      mat_test = strfind(d(i).name, '.mat');
      chi_test = strfind(d(i).name, 'chi');
      ic_test  = strfind(d(i).name, 'ic.');

      if ~isempty(mat_test) & ~isempty(chi_test)

         ID = d(i).name(1:mat_test-1);
         if isChipod
             if ID(end) == '1' | strcmpi(ID(end-3:end), '1_ic') % sensor T1
                 sensor = 1;
             end
             if ID(end) == '2' | strcmpi(ID(end-3:end), '2_ic') % sensor T1
                 sensor = 2;
             end
         else
             sensor = 1; % gusT has only one sensor
         end

         if ~isempty(strfind(ID, 'p'))
             isPitotEstimate = 1;
         else
             isPitotEstimate = 0;
         end

         disp(' ');
         disp(['----------> adding ' ID ]);
         load([dirname ID '.mat'])

         % find desired time range
         iiTrange = find( chi.time >= time_range(1) & chi.time<= time_range(2) );

         % extract valid time range here so that histograms work
         ff = fieldnames(chi);
         for nn = 1:length(ff)
             try
                 chi.(ff{nn}) = chi.(ff{nn})(iiTrange);
             catch ME;
             end
         end

         try
             load ../proc/T_m.mat
             Smean = interp1(T1.time, (T1.S + T2.S)/2, chi.time);
             chi.S = Smean;
         catch ME
             Smean = 35*ones(size(chi.time));
         end
         rho = sw_pden(Smean, chi.T, ChipodDepth, 0);
         cp = sw_cp(Smean, chi.T, ChipodDepth);

         %___________________NaN out after sensor death______________________
         if sensor == 1 % sensor T1 or gusT T
             death = find(chi.time > T1death, 1, 'first');
             if ~isempty(death)
                 disp(['NaNing out sensor T1 after it died on ' datestr(T1death)])
                 chi.chi(death:end) = NaN;
                 chi.eps(death:end) = NaN;
             end
         end

         if isChipod
             if sensor == 2 % sensor T2
                 death = find(chi.time > T2death, 1, 'first');
                 if ~isempty(death)
                     disp(['NaNing out sensor T2 after it died on ' datestr(T2death)])
                     chi.chi(death:end) = NaN;
                     chi.eps(death:end) = NaN;
                 end
             end
         end

         %____________________NaN out specific time ranges as necessary____________
         % for temp sensor
         if ~isempty(nantimes{sensor})
             for tt = 1:size(nantimes{sensor}, 1)
                 chi.chi(find_approx(chi.time, nantimes{sensor}(tt, 1), 1): ...
                         find_approx(chi.time, nantimes{sensor}(tt, 2), 1)) = NaN;
             end
             chi.eps(isnan(chi.chi)) = NaN;
         end

         if isPitotEstimate & ~isempty(nantimes{3})
              for tt = 1:size(nantimes{3}, 1)
                 chi.chi(find_approx(chi.time, nantimes{3}(tt, 1), 1): ...
                         find_approx(chi.time, nantimes{3}(tt, 2), 1)) = NaN;
             end
             chi.eps(isnan(chi.chi)) = NaN;
         end

         chi.Kt = 0.5 * chi.chi ./ chi.dTdz.^2;
         chi.Jq = -rho .* cp .* chi.Kt .* chi.dTdz;

         if do_plot
             hfig = CreateFigure;
             Histograms(chi, hfig, normstr, 'raw');
         end

         if do_mask
             % speed mask could change depending on estimate
             if strcmpi(mask_spd_initial, '')
                 mask_spd = ID(5);
             end

             if mask_spd == 'm' & ~exist('vel_m', 'var')
                 load ../input/vel_m.mat
                 vel = vel_m;
                 disp('masking using mooring speed');
             elseif mask_spd == 'p' & ~exist('vel_p', 'var')
                 load ../input/vel_p.mat
                 vel = vel_p;
                 disp('masking using pitot speed');
             end
             spdmask = interp1(vel.time, vel.spd, chi.time);

             % histograms for speed-based masking
             if do_plot & exist('../proc/temp.mat', 'file')
                 load ../proc/temp.mat
                 load ../input/vel_m.mat
                 CreateFigure;

                 dt = diff(T.time(1:2))*86400;
                 trange = [find_approx(T.time, time_range(1)):find_approx(T.time, time_range(2))];
                 subplot(211)
                 histogram(abs(T.a_vel_x(trange)), 'normalization' ,'pdf', ...
                           'displaystyle', 'stairs', 'linewidth', 1.5);
                 hold on;
                 histogram(hypot(T.a_vel_x(trange), T.a_vel_y(trange)), ...
                           'normalization' ,'pdf', 'displaystyle', 'stairs', 'linewidth', 1.5);
                 histogram(abs(T.a_vel_z(trange)), 'normalization' ,'pdf', ...
                           'displaystyle', 'stairs', 'linewidth', 1.5);

                 trange2 = [find_approx(vel_m.time, time_range(1)):find_approx(vel_m.time, time_range(2))];
                 histogram(vel_m.spd(trange2), 'normalization' ,'pdf', ...
                           'displaystyle', 'stairs', 'linewidth', 1.5);
                 plot([1 1]*min_spd, ylim, 'k--');
                 xlim([-0.25 1]*0.7*max(vel_m.spd))
                 set(gca, 'XTick', sort([min_spd get(gca, 'XTick')]));
                 xlabel('m/s')
                 ylabel('PDF')
                 legend('|v_x|', '|v_x + v_y|', '|v_z|', 'background flow', 'min. spd criterion');
                 dt2 = diff(vel_m.time(1:2))*86400/60;
                 title([ID(5:end) ' | ' num2str(round(dt)) 's avg for v_* | ' ...
                        num2str(round(dt2)) 'minutes \Deltat for background'])

                 subplot(212)
                 histogram(log10(chi.eps), 'normalization', 'pdf', 'displaystyle', 'stairs', 'linewidth', 1.5);
                 hold on;
                 histogram(log10(chi.eps(spdmask < min_spd)), ...
                           'normalization', 'pdf', 'displaystyle', 'stairs', 'linewidth', 1.5);
                 ylabel('PDF')
                 xlabel('log_{10} \epsilon')
                 xlim([-14 5])
                 legend('raw', ['background flow < ' num2str(min_spd) 'm/s']);

                 print(gcf,['../pics/velocity-masking-' ID(5:end) '.png'],'-dpng','-r200','-painters')
                 savefig(gcf,['../pics/velocity-masking-' ID(5:end) '.fig'])
             end

             chi = ApplyMask(chi, abs(chi.dTdz), '<', min_dTdz, 'Tz');
             if do_plot, Histograms(chi, hfig, normstr, 'Tz'); end

             chi = ApplyMask(chi, chi.N2, '<', min_N2, 'N2');
             if do_plot, Histograms(chi, hfig, normstr, 'N2'); end

             chi = ApplyMask(chi, chi.spd, '<', min_inst_spd, 'inst speed');
             chi = ApplyMask(chi, spdmask, '<', min_spd, 'background flow');

             % additional Tz masking?
             if ~isempty(Tz)
                 if mask_dTdz == 'i'
                     % choose appropriate internal stratification for sensor
                     Tz.Tz = Tz.(['Tz' ID(7)']);
                 end

                 Tzmask = interp1(Tz.time, Tz.Tz, chi.time);
                 chi = ApplyMask(chi, abs(Tzmask), '<', 1e-4, ...
                                 ['Additional Tz_' mask_dTdz]);
             end

             if do_plot
                 Histograms(chi, hfig, normstr, 'all masks');

                 figure(hfig)
                 subplot(221); legend(gca, 'show'); title(ID(5:end));
                 subplot(222); legend(gca, 'show'); title(ID(5:end));
                 subplot(223); legend(gca, 'show'); title(ID(5:end));
                 subplot(224); legend(gca, 'show'); title(ID(5:end));

                 print(gcf,['../pics/histograms-masking-' ID '.png'],'-dpng','-r200','-painters')
                 savefig(gcf,['../pics/histograms-masking-' ID '.fig'])
             end
         end

         % convert averaging window from seconds to points
         ww =  round(avgwindow/(diff(chi.time(1:2))*3600*24));

         if isempty(ic_test)
             % deglitch chi and eps before
             % calculating Jq and Kt
             % not required for IC estimate because that is already
             % an averaged estimate
             disp('Deglitch... itch... tch... ch')
             tic;
             chi.chi = deglitch(chi.chi, ww, 2,'b');
             chi.eps = deglitch(chi.eps, ww, 2, 'b');
             toc;

             % get list of all fields to average
             ff = fields(chi);

             %% average data
             disp('Running moving average')
             tic;
             for f = 1:length(ff)  % run through all fields in chi
                 if ( length(chi.(ff{f})) == length(chi.time) )
                     Turb.(ID).(ff{f}) = moving_average( chi.(ff{f}), ww, ww );
                 end
             end
             toc;
         else
             Turb.(ID) = chi;
         end

         if do_plot
             if ~exist('hfig2', 'var'), hfig2 = CreateFigure; end
             Histograms(Turb.(ID), hfig2, 'pdf', ID);
         end

      end
   end

   if do_plot
       figure(hfig2)
       subplot(221); legend(gca, 'show'); title(['Final ' num2str(avgwindow/60) ' min mean']);
       subplot(222); title(['Final ' num2str(avgwindow/60) ' min mean']);

       print(gcf,['../pics/histograms-final.png'],'-dpng','-r200','-painters')
       savefig(gcf,['../pics/histograms-final.fig'])
   end

   Turb.do_mask = do_mask;
   Turb.mask_dTdz = mask_dTdz;
   Turb.min_dTdz = min_dTdz;
   Turb.min_spd = min_spd;
   Turb.avgwindow = avgwindow;
   Turb.min_inst_spd = min_inst_spd;
   Turb.min_N2 = min_N2;
   Turb.mask_spd = mask_spd_initial;

   %---------------------add readme----------------------
   Turb.readme = {...
         '------------------sub-fields--------------------'; ...
         'time    :  matlab time'; ...
         'chi     :  chi (10 min averaged deglitched)'; ...
         'eps     :  dissipation (10 min averaged deglitched)'; ...
         'dTdz    :  vertical temperature gradient used for processing'; ...
         'N2      :  buoyancy frequency used for processing'; ...
         'spd     :  flow pass sensor (for pumped moorings contains ACC-information)'; ...
         'T       :  temperature used for calculation'; ...
         'S       :  salinity used for calculation (default 35)'; ...
         'mask    :  values that got flagged during processing due to minTz or speed'; ...
         ''; ...
         '------------------main-fields--------------------'; ...
         '..._ic  :  inertial convective range estimate (no subscript viscous convective range)'; ...
         'mm1     :  mooring speed  - mooring N2              - temp_sensor 1 '; ...
         'mm2     :  mooring speed  - mooring N2              - temp_sensor 2 '; ...
         'mi11    :  mooring speed  - internal N2 (T1 sensor) - temp_sensor 1 '; ...
         'mi112   :  mooring speed  - internal N2 (combination T1 and T2 sensor) - temp_sensor 1 '; ...
         'mi22    :  mooring speed  - internal N2 (T2 sensor) - temp_sensor 2 '; ...
         'pm1     :  Pitot speed    - mooring N2              - temp_sensor 1 '; ...
         'pm2     :  Pitot speed    - mooring N2              - temp_sensor 2 '; ...
         'pi11    :  Pitot speed    - internal N2 (T1 sensor) - temp_sensor 1 ' ...
         };

%_____________________save combined structure______________________
   save([savedir '/Turb.mat'], 'Turb');
   if length(runname) ~= 0
       system(['ln -s ' savedir '/Turb.mat ' savedir ...
               '../combined/' runname '.mat']);
   end
end

%_____________________comparison plot______________________
if do_plot
   
   load([basedir '/proc/' runname '/Turb.mat']);
   ff = fields(Turb);
   ff = {ff{1:end-1}}'; % remove readme structure

   fig = CreateFigure;
    % fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
    %         'Papersize',[30 20],'PaperPosition',[0 0 30 20])
    
         [ax, ~] = create_axes(fig, 4, 1, 0);
      
         col = get(groot,'DefaultAxesColorOrder');
         col = cat(1, col, col./2, (1-col)/2+col); % colorscale extension
         
         a=1;
         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})), continue; end
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).chi, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['\chi [K^2/s]'], 1);
         set(ax(a), 'Yscale', 'log');
         yl = [1e-10 1e-2];
         ylim(ax(a), yl);
         
         a=2;
         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})) == 1, continue; end
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).eps, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['\epsilon [m^2/s^3]'], 1);
         set(ax(a), 'Yscale', 'log');
         yl = [1e-10 1e-2];
         ylim(ax(a), yl);
         
         
         a=3;
         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})) == 1, continue; end
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).N2, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['N^2 [s^{-2}]'], 1);
         set(ax(a), 'Yscale', 'log');
         yl = [1e-6 1e-2];
         ylim(ax(a), yl);
         

         a=4;
         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})) == 1, continue; end
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).spd, 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         t = text_corner(ax(a), ['|u| [m/s]'], 1);
         datetick(ax(a), 'keeplimits');
         yl = [0 2];
         ylim(ax(a), yl);

         linkaxes(ax, 'x');
         
            
   %---------------------histogram plots----------------------
      squeeze_axes(ax, .8, 1)

      pos = get(ax(1), 'Position');
      axh(1) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.01 pos(4)] );
      
      pos = get(ax(2), 'Position');
      axh(2) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.01 pos(4)] );

         a=1;
         hold(axh(a),'all');
         yl = log10(get(ax(a), 'Ylim'));
         bins = yl(1):diff(yl)/100:yl(2);
         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})) == 1, continue; end
            [Nchi,~] = histcounts( log10(Turb.(ff{f}).chi) , bins);
            pj = f; p(pj) = plot(axh(a), Nchi , bins(1:end-1)+diff(bins(1:2)*.5), 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         ylim(axh(a), yl);
         set(axh, 'Yticklabel', {}, 'Xticklabel', {})
      
         a=2;
         hold(axh(a),'all');
         yl = log10(get(ax(a), 'Ylim'));
         bins = yl(1):diff(yl)/100:yl(2);
         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})) == 1, continue; end
            [Nchi,~] = histcounts( log10(Turb.(ff{f}).eps) , bins);
            pj = f; p(pj) = plot(axh(a), Nchi , bins(1:end-1)+diff(bins(1:2)*.5), 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         ylim(axh(a), yl);
         set(axh, 'Yticklabel', {}, 'Xticklabel', {})
         

   %---------------------legend----------------------
      pos = get(ax(3), 'Position');
      axl = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.01 pos(4)] );
      hold(axl,'on');
      

         for f = 1:length(ff)
            if ~isstruct(Turb.(ff{f})) == 1, continue; end
            pj = f; p(pj) = plot(axl, [0 1] ,[0 1], 'color', [col(pj,:) 1], 'Linewidth', 1);
         end
         legend(p, ff);
         set(axl, 'visible', 'off')
         ylim(axl,[-1 -.5])


   %---------------------save imagage----------------------

   print(gcf,'../pics/Compare_Turb.png','-dpng','-r200','-painters')
   savefig(gcf,'../pics/Compare_Turb.fig')
   
end


% Examples of using TestMask and DebugPlots to check masking
% TestMask(chi, abs(chi.dTdz), '<', [1e-4, 3e-4, 1e-3], 'Tz');
% t0 = datenum(2014, 01, 01);
% t1 = datenum(2014, 03, 01);
% DebugPlots([], t0, t1, chi, 'raw', 1)
% chi1 = ApplyMask(chi, abs(chi.dTdz), '<', 3e-4, 'T_z');
% DebugPlots([], t0, t1, chi1, 'T_z > 3e-4', 1)
% load ../proc/temp.mat
% load ../proc/Turb.mat
