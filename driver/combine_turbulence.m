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
   min_N2 = 1e-6;
   min_dTdz = 1e-4;
   min_spd = 0.05;
   min_inst_spd = min_spd; % min instantaneous speed past sensor
   mask_dTdz = 'i'; % '' for none, 'm' for mooring, 'i' for internal
                    % in addition to whats in chi.dTdz
   mask_inst_spd = 1; % estimates are crappy if sensor isn't moving
                      % enough.
                      % screws the spectrum calculation...
   mask_spd = ''; % masking such that background flow flushes
                  % sensed water volume
                  % 'm' for mooring, 'p' for pitot,
                  % '' to choose based on what was used in chi estimate
   avgwindow = 600; % averaging window in seconds


   % if you want to restrict the time range that should be combined
   % use the following
   % This restricts time range for BOTH sensors on chipods
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

   % if on chipod one sensor dies earlier, specify that time here.
   T1death = datenum(2060, 1, 1, 0, 0, 0);
   T2death = datenum(2060, 1, 1, 0, 0, 0);


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
   d = dir([basedir '/proc/chi/']);
   % eliminate directories
   d = d(~[d(:).isdir]);

   % run through all chi files
   for i = 1:length(d)

      mat_test = strfind(d(i).name, '.mat');
      chi_test = strfind(d(i).name, 'chi');
      ic_test  = strfind(d(i).name, 'ic.');

      if ~isempty(mat_test) & ~isempty(chi_test)

         ID = d(i).name(1:mat_test-1);
         disp(['----------> adding ' ID ]);
         load([basedir '/proc/chi/' ID '.mat'])
         
         % find desired time range
         iiTrange = find( chi.time>= time_range(1) & chi.time<= time_range(2) );

         ChipodDepth = round(nanmean(chi.depth(iiTrange))); % in metres

         % read in mooring salinity & calculate rho, cp
         if ~exist('Smean', 'var')
             try
                 load ../proc/T_m.mat
                 Smean = (T1.S .* abs(nanmean(T2.z)-ChipodDepth) ...
                          + T2.S .* abs(nanmean(T1.z)-ChipodDepth)) ...
                         ./ abs(nanmean(T1.z)-nanmean(T2.z));
                 Smean = interp1(T1.time, Smean, chi.time);
                 chi.S = Smean;
             end
             rho = sw_pden(Smean, chi.T, ChipodDepth, 0);
             cp = sw_cp(Smean, chi.T, ChipodDepth);
         end

         if do_mask
             percent_mask_Tz = sum(abs(chi.dTdz(iiTrange)) < min_dTdz)/length(chi.time(iiTrange))*100;
             percent_mask_N2 = sum(chi.N2(iiTrange) < min_N2)/length(chi.time(iiTrange))*100;
             percent_mask_inst_spd = sum(chi.spd(iiTrange) < min_inst_spd)/length(chi.time(iiTrange))*100;
             disp([' Tz will mask ' num2str(percent_mask_Tz, '%.2f') ...
                   ' % of estimates'])
             disp([' N2 will mask ' num2str(percent_mask_N2, '%.2f') ...
                   ' % of estimates'])
             disp([' Inst speed will mask ' num2str(percent_mask_inst_spd, '%.2f') ...
                   ' % of estimates'])

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
             percent_mask_spd = sum(spdmask(iiTrange) < min_spd)/length(chi.time(iiTrange))*100;
             disp([' speed will mask ' num2str(percent_mask_spd, '%.2f') ...
                   '% of estimates'])

             full_mask = (abs(chi.dTdz) < min_dTdz) ...
                 | (chi.N2 < min_N2) ...
                 | (chi.spd < min_inst_spd) ...
                 | (spdmask < min_spd);

             if ~isempty(Tz) % additional Tz masking?
                 if mask_dTdz == 'i'
                     % choose appropriate internal stratification for sensor
                     Tz.Tz = Tz.(['Tz' ID(7)']);
                 end

                 Tzmask = interp1(Tz.time, Tz.Tz, chi.time);
                 percent_mask_dTdz = sum(abs(Tzmask(iiTrange)) < min_dTdz)/ ...
                     length(chi.time(iiTrange))*100;
                 disp([' dTdz_' mask_dTdz ' will mask ' num2str(percent_mask_dTdz, '%.2f') ...
                       ' % of estimates'])

                 full_mask = full_mask | (abs(Tzmask) < min_dTdz);
             end
         end

         % NaN out some chi estimates based on min_dTz, min_spd
         if do_mask
             chi.chi(full_mask) = NaN;
             chi.eps(full_mask) = NaN;
             chi.mask = chi.mask | ~full_mask;
         end

         % convert averaging window from seconds to points
         ww =  round(avgwindow/(diff(chi.time(1:2))*3600*24));

         if isempty(ic_test)
             % deglitch chi and eps before
             % calculating Jq and Kt
             % not required for IC estimate because that is already
             % an averaged estimate
             chi.chi = deglitch(chi.chi, ww, 2,'b');
             chi.eps = deglitch(chi.eps, ww, 2, 'b');
         end

         if isChipod
             if ID(end) == '1' % sensor T1
                 death = find(chi.time > T1death, 1, 'first');
                 if ~isempty(death)
                     disp(['NaNing out sensor T1 after it died on ' datestr(T1death)])
                     chi.chi(death:end) = NaN;
                     chi.eps(death:end) = NaN;
                 end
             end

             if ID(end) == '2' % sensor T2
                 death = find(chi.time > T2death, 1, 'first');
                 if ~isempty(death)
                     disp(['NaNing out sensor T2 after it died on ' datestr(T2death)])
                     chi.chi(death:end) = NaN;
                     chi.eps(death:end) = NaN;
                 end
             end
         end

         chi.Kt = 0.5 * chi.chi ./ chi.dTdz.^2;
         chi.Jq = -rho .* cp .* chi.Kt .* chi.dTdz;

         % get list of all fields to average
         ff = fields(chi);

         %% average data
         disp('Running moving average')
         for f = 1:length(ff)  % run through all fields in chi
             if ( length(chi.(ff{f})) == length(chi.time) )
                 Turb.(ID).(ff{f}) = moving_average( chi.(ff{f})(iiTrange), ww, ww );
             end
         end
      end
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
   save([basedir '/proc/Turb.mat'], 'Turb');

    
end

%_____________________comparison plot______________________
if do_plot
   
   load([basedir '/proc/Turb.mat']);
   ff = fields(Turb);
   ff = {ff{1:end-1}}'; % remove readme structure

    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 20],'PaperPosition',[0 0 30 20])
    
         [ax, ~] = create_axes(fig, 4, 1, 0);
      
         col = get(groot,'DefaultAxesColorOrder');
         col = cat(1, col, col./2); % colorscale extension
         
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
   
end
