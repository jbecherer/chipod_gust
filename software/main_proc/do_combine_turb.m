function [] = do_combine_turb( basedir, savedir, CP, do_combine, do_plot, do_mask, save_fig, is_visible) 
%%    [] = do_combine_turb( basedir, savedir, CP, do_combine, do_plot, do_mask, save_fig, is_visible) 
%
%        All the fuction of combine turbulence in one routine
%
%        INPUT
%           basedir     :  unit directory
%           savedir     :  where shall the structure be saved
%           CP          :  flagging parameters generated with default_parameters_combine_turbulence.m
%           do_combine  :  1: do_combining 0:don't
%           do_plot     :  generates plots (1 yes 0 no)
%           do_mask     :  apply mask
%           save_fig    :  save fig files
%           is_visible  :  visibilty of figures 'on' or 'off'
%
%
%   created: 
%        Thu Dec 28 18:18:54 PST 2017

%_____________________defaults______________________
if nargin < 8
   is_visible = 'on';
end
if nargin < 7
   save_fig   = 0;
end
if nargin < 6
   do_mask = 1;
end
if nargin < 5
   do_plot  =  1;
end
if nargin < 4
   do_combine;
end
if nargin < 3
   [CP] = default_parameters_combine_turbulence(basedir);
end
if nargin < 2
   savedir = [basdir 'proc/'];
end

wdafile = [basedir 'input' filesep 'dTdz_w.mat'];
motionfile = [basedir 'proc' filesep 'motion.mat'];

% determine if chipod or gusT
load([basedir '/calib/header.mat'])
if isfield(head.coef, 'T1')
    CP.isChipod = 1;
else
    CP.isChipod = 0;
end

if do_mask, Tz = determine_additional_Tz_mask(CP); end

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

         ID = d(i).name(5:mat_test-1);
         if ~CP.pflag.proc.(ID)
             disp([ ID ' is disabled!']);
             continue;
         end

         disp(' ');
         disp(['----------> adding ' ID ]);
         clear chi;
         load([dirname 'chi_' ID '.mat'])

         fds = fields(chi);
         for f = 1:length(fds)
            if size(chi.(fds{f}),1) ~= 1
               chi.(fds{f}) = chi.(fds{f})';
            end
         end

         CP = process_estimate_ID(CP, ID);

         % do winters dasaro estimate?
         do_wda = (exist(wdafile, 'file') || isfield(chi, 'wda')) ...
                  && ~contains(ID, 'ic') && CP.pflag.master.winters_dasaro;

         % convert averaging window from seconds to points
         ww = round(CP.avgwindow/(diff(chi.time(1:2))*3600*24));
         dw = round(CP.deglitch_window/(diff(chi.time(1:2))*3600*24));

         chi = truncate_time(chi, CP.time_range);

         try
             load([basedir '/proc/T_m.mat']);
             chi.S = interp1(T1.time, (T1.S + T2.S)/2, chi.time);
             compare_min_dTdz_N2_against_SBE_specs(nanmedian(abs(T1.z - T2.z)), ...
                                                   CP, ID);
         catch ME
             chi.S = 35*ones(size(chi.time));
         end

         % add nans when sensors die or glitch
         chi = add_nans(CP, chi);

         chi.Kt = 0.5 * chi.chi ./ chi.dTdz.^2;
         chi.Jq = -1025 .* 4200 .* chi.Kt .* chi.dTdz;

         % save unmasked chi structure for later use
         chiold = chi;

         if do_plot
             hfig = CreateFigure(is_visible);
             hfig.Name = ['Histograms: effect of masking for ' ID];
             Histograms(chi, hfig, CP.normstr, ID, 'raw');

             if ~exist('hfraw', 'var')
                 hfraw = CreateFigure(is_visible);
                 hfraw.Name = ['Histograms: compare different 1 sec estimates'];
             end
             Histograms(chi, hfraw, 'pdf', ID, ID);

             if ~exist('hfstrat', 'var')
                 hfstrat = CreateFigure(is_visible);
                 hfstrat.Name = ['Histograms: background stratification for ' ID];
                 shown_Tz = '';
             end

             % If we haven't shown histograms of the current Tz estimate, do so.
             if isempty(strfind(shown_Tz, ID(2)))
                 StratHist(hfstrat, chi, ID);
                 set(0, 'currentfigure', hfstrat); subplot(122);
                 hplt = plot(CP.avgwindow/60*[1, 1], ylim, 'k--');
                 legend(hplt, 'averaging window')
                 shown_Tz = [shown_Tz ID(2)];
             end
         end

         if do_mask
             [spdmask, addspdmask] = determine_speed_masks(basedir, ID, CP, chi);

             if CP.mask_flushing

                 motion = process_motion_file(motionfile);

                 if ~exist('badMotion', 'var') | length(badMotion)~=length(chi.chi)
                     badMotion = make_flushing_mask(motion, mask_spd, vel, do_plot);
                     if do_plot
                         print(gcf, [basedir 'pics' filesep 'angles.png'], '-r200', '-painters', '-bestfit')
                         if save_fig, set(gcf, 'visible', 'on'); savefig(gcf, [basedir 'pics' filesep 'angles.fig']); end
                     end
                 end

                 chi = ApplyMask(chi, badMotion, '=', 1, 'volume not being flushed');
                 if do_plot, Histograms(chi, hfig, CP.normstr, (ID), 'volume flushed'); end
             end

             [chi, percentage] = ApplyMask(chi, chi.N2, '<', CP.min_N2, 'N2');
             chi.stats.N2_mask_percentage = percentage;
             if percentage > 0.5
                 perlabel = [' -' num2str(percentage, '%.1f') '%'];
                 if do_plot, Histograms(chi, hfig, CP.normstr,(ID), ['N2' perlabel]); end
             end

             [chi, percentage] = ApplyMask(chi, chi.spd, '<', CP.min_inst_spd, 'inst speed');
             chi.stats.inst_speed_mask_percentage = percentage;
             [chi, percentage] = ApplyMask(chi, spdmask, '<', CP.min_spd, 'background flow');
             chi.stats.background_flow_mask_percentage = percentage;

             if CP.additional_mask_spd ~= ''
                 [chi, percentage] = ApplyMask(chi, addspdmask, '<', CP.min_spd, 'background flow');
                 chi.stats.additional_background_flow_mask_percentage = percentage;
             end

             % remove values greater than thresholds
             [chi, chi.stats.max_chi_percentage] = ApplyMask(chi, chi.chi, '>', CP.max_chi, 'max_chi');
             [chi, chi.stats.max_eps_percentage] = ApplyMask(chi, chi.eps, '>', CP.max_eps, 'max_eps');

             if isempty(ic_test)
                 % deglitch chi and eps before calculating Jq and Kt
                 % not required for IC estimate because that is already
                 % an averaged estimate
                 disp('Deglitch... itch... tch... ch')
                 tic;
                 chi.chi = 10.^deglitch(log10(chi.chi), dw, CP.deglitch_nstd, 'b');
                 chi.eps(isnan(chi.chi)) = nan;
                 chi.eps = 10.^deglitch(log10(chi.eps), dw, CP.deglitch_nstd, 'b');
                 chi.chi(isnan(chi.eps)) = nan;
                 toc;
             end

             if isempty(ic_test) & isfield(chi, 'spec_area') & isfield(chi, 'spec_floor')
                 % account for two possible versions of the code
                 if ~isfield(chi, 'time_floor')
                     spec_floor = nanmedian(chi.spec_floor);
                 else
                     if length(chi.time_floor)>1
                        spec_floor = interp1(chi.time_floor, chi.spec_floor, chi.time);
                     else
                        spec_floor = chi.spec_floor;
                     end
                 end

                 if ~exist('hspecfig', 'var')
                     hspecfig = CreateFigure(is_visible);
                     hspecfig.Position(3) = 2000;
                     hspecfig.Position(4) = 1065;
                     hspec1 = subplot(2,2,[1,2]); hold on;
                     hspec2 = subplot(223); hold on;
                     hspec3 = subplot(224); hold on;

                     shown_sensors_noise_floor = [];
                 end

                 spec_floor_mask = chi.spec_area < CP.factor_spec_floor * spec_floor * nanmean(chi.nfft);

                 do_spec_area = ~any(shown_sensors_noise_floor == CP.sensor);
                 make_noise_floor_histograms(ID, chi, hspec1, hspec2, hspec3, ...
                                             spec_floor, spec_floor_mask, do_spec_area)
                 if do_spec_area, shown_sensors_noise_floor = [shown_sensors_noise_floor, CP.sensor]; end

                 [chi, chi.stats.spec_floor_percentage] = ApplyMask(chi, ...
                                                                   spec_floor_mask & isnan(chi.chi), '=', 1, ...
                                                                   'spec_floor', [], CP.noise_floor_fill_value);
             end

             % filter out bad fits using fitting statistics
             statsname = [basedir filesep 'proc' filesep 'chi' filesep 'stats' ...
                         filesep 'chi_' ID '_stats.mat'];
             if ~contains(ID, '_ic')
                 if exist(statsname, 'file')
                     load(statsname);
                     stats = truncate_time(stats, CP.time_range);
                     assert(isequal(stats.time, chi.time));

                     [chi, chi.stats.nfreq_percentage] = ...
                         ApplyMask(chi, stats.n_freq, '<', CP.min_n_freq, 'min n_freq');

                     if CP.mask_ic_fits
                         [chi, chi.stats.ic_fit_percentage] = ...
                             ApplyMask(chi, stats.k_stop < stats.ki, '=', 1, ...
                                       'fully IC fit');
                         perlabel = [' -' num2str(chi.stats.nfreq_percentage + chi.stats.ic_fit_percentage, '%.1f') '%'];
                         if do_plot, Histograms(chi, hfig, CP.normstr, (ID), ['bad fits' perlabel]); end
                     end

                 else
                     disp(['Combined stats file does not exist. Cannot filter out bad fits.'])
                 end
             end
             % obtain Kt, Jq using Winters & D'Asaro methodology
             if do_wda
                 ticwda = tic;
                 disp('Processing Winters & D''Asaro estimate. Takes 120s for 1 year.');
                 if ~exist(wdafile, 'file') && isfield(chi, 'wda')
                     % backward compatibility
                     disp(['    ' wdafile ' not found. \n' ...
                           '     Using chi.wda to do Winters & D''Asaro estimate'])
                     chi.wda = process_wda_estimate(chi, chi.wda);
                 else
                     if ~exist('Tz_w', 'var') & exist(wdafile, 'file')
                         % only needs to be loaded once
                         load(wdafile);
                     end
                     if CP.isChipod
                         chi.wda = process_wda_estimate(chi, Tz_w.(['wda' num2str(CP.sensor)]));
                     else
                         chi.wda = process_wda_estimate(chi, Tz_w.wda);
                     end
                 end
                 toc(ticwda);

                 if ~exist(wdafile, 'file') && isfield(chi, 'wda')
                     if CP.wda_Tz_sign == 'm'
                         % determine sign of /vertical/ heat flux using mooring gradient
                         fname = [basedir 'input' filesep 'dTdz_m.mat'];
                         if exist(fname, 'file')
                             Tz = load(fname);
                             Tz = Tz.Tz_m;
                         else
                             error([fname 'does not exist. Specify CP.wda_Tz_sign ' ...
                                    'properly.'])
                         end
                     elseif CP.wda_Tz_sign == 'i'
                         fname = [basedir 'input' filesep 'dTdz_i.mat'];
                         if exist(fname, 'file')
                             Tz = load(fname);
                             Tz = Tz.Tz_i;
                             Tz.Tz = Tz.(['Tz' num2str(CP.sensor)]);
                         else
                             error([fname 'does not exist. Specify CP.wda_Tz_sign ' ...
                                    'properly.'])
                         end
                     end

                     sgn = get_wda_sign(chi.wda.time, Tz, CP.min_dTdz);
                     if CP.wda_Tz_sign == 'm', chi.wda.sign_used = 'mooring'; end
                     if CP.wda_Tz_sign == 'i', chi.wda.sign_used = 'internal'; end
                 else
                     if CP.wda_Tz_sign == 'm'
                         sgn = Tz_w.sgn_moor;
                         chi.wda.sign_used = 'mooring';
                     elseif CP.wda_Tz_sign == 'i'
                         sgn = Tz_w.(['sgn_int' num2str(CP.sensor)]);
                         chi.wda.sign_used = 'internal';
                     end
                 end

                 chi.wda.sgn = sgn;
                 chi.wda.Jq = -abs(chi.wda.Jq) .* chi.wda.sgn;
                 chi.wda.dTdz = abs(chi.wda.dTdz) .* chi.wda.sgn;

                 % add in molecular diffusivity
                 chi.wda.Kt = chi.wda.Kt + sw_tdif(interp1(chi.time, chi.S, chi.wda.time), ...
                                                   interp1(chi.time, chi.T, chi.wda.time), ...
                                                   CP.depth);

                 % add sorted gradient to stratification histograms
                 if do_plot && isempty(strfind(shown_Tz, 'w'))
                     shown_Tz = [shown_Tz 'w'];
                     StratHist(hfstrat, chi.wda, 'ww');
                 end
             end

             % dT/dz has to happen after I use chi to get Winters & D'Asaro
             % estimates of Kt, Jq!
             [chi, percentage] = ApplyMask(chi, abs(chi.dTdz), '<', CP.min_dTdz, 'Tz');
             chi.stats.dTdz_mask_percentage = percentage;
             perlabel = [' -' num2str(percentage, '%.1f') '%'];
             if do_plot, Histograms(chi, hfig, CP.normstr, (ID), ['|Tz| > ' num2str(CP.min_dTdz, '%.1e') perlabel]); end

             % additional Tz masking?
             if CP.additional_mask_dTdz ~= ''
                 if CP.additional_mask_dTdz == 'i'
                     % choose appropriate internal stratification for sensor
                     Tz.Tz = Tz.(['Tz' ID(3)']);
                 end

                 Tzmask = interp1(Tz.time, Tz.Tz, chi.time);
                 [chi, percentage] = ApplyMask(chi, abs(Tzmask), '<', 1e-4, ...
                                               ['Additional Tz_' CP.additional_mask_dTdz]);
                 chi.stats.additional_dTdz_mask_percentage = percentage;
                 perlabel = [' -' num2str(percentage, '%.1f') '%'];
                 if do_plot, Histograms(chi, hfig, CP.normstr, (ID), ...
                                        ['Additional Tz_' CP.additional_mask_dTdz perlabel]); end
             end

             if do_plot
                 set(0, 'currentfigure', hfig);
                 set(hfig, 'DefaultLegendBox', 'off');
                 subplot(221); legend(gca, 'show'); title((ID));
                 subplot(222); legend(gca, 'show'); title((ID));
                 subplot(223); legend(gca, 'show'); title((ID));
                 subplot(224); legend(gca, 'show'); title((ID));

                 print(gcf,[basedir '/pics/histograms-masking-' ID '.png'],'-dpng','-r200','-painters')
                 if save_fig,set(gcf, 'visible', 'on'); savefig(gcf,[basedir '/pics/histograms-masking-' ID '.fig']); end
             end
         end

         if isempty(ic_test)
             % get list of all fields to average
             ff = fields(chi);

             %% average data
             disp('Running moving average')
             tic;
             for f = 1:length(ff)  % run through all fields in chi
                 if ( length(chi.(ff{f})) == length(chi.time) )
                     if strcmpi(ff{f}, 'Kt') | strcmpi(ff{f}, 'Jq'), continue; end
                     Turb.(ID).(ff{f}) = moving_average( chi.(ff{f}), ww, ww , CP.avgvalid);
                 else
                     Turb.(ID).(ff{f}) = chi.(ff{f});
                 end
             end
             toc;
         else
             Turb.(ID) = chi;
         end

         % reapply dTdz filter
         % (because we divide by this dTdz for Kt, Jq)
         Turb.(ID) = ApplyMask(Turb.(ID), abs(Turb.(ID).dTdz), '<', CP.min_dTdz, 'avg dTdz');

         % Tz (min_dTdz)-crossing filter
         Tz_min_cross = generate_min_dTdz_crossing_mask(Turb.(ID).dTdz, ...
                                                        CP.min_dTdz, 0);
         Turb.(ID) = ApplyMask(Turb.(ID), Tz_min_cross, '=', 1, 'min_Tz crossing');

         % recalculate using averaged quantities
         % if we average over a time period greater than
         % sampling period of dTdz, this estimate will differ!
         Turb.(ID).Kt = 0.5 * Turb.(ID).chi ./ Turb.(ID).dTdz.^2 + ...
             sw_tdif(interp1(chi.time, chi.S, Turb.(ID).time), ...
                     Turb.(ID).T, CP.depth);
         Turb.(ID).Jq = -1025 .* 4200 .* Turb.(ID).Kt .* Turb.(ID).dTdz;

         [Turb.(ID), Turb.(ID).stats.max_Kt_percentage] = ApplyMask(Turb.(ID), Turb.(ID).Kt, '>', CP.max_Kt, 'max_Kt');
         [Turb.(ID), Turb.(ID).stats.max_Jq_percentage] = ApplyMask(Turb.(ID), abs(Turb.(ID).Jq), '>', CP.max_Jq, 'max_Jq');
            
         if do_wda
             % get list of all fields to average
             ff = fields(chi.wda);

             %% average data
             % disp('WDA: Running moving average')
             % wwwda = round(CP.avgwindow/(diff(chi.wda.time(1:2))*24*3600));
             % tic;
             % for f = 1:length(ff)  % run through all fields in chi
             %     if ( length(chi.wda.(ff{f})) == length(chi.wda.time) )
             %         Turb.(ID).wda.(ff{f}) = moving_average( chi.wda.(ff{f}), wwwda, wwwda , CP.avgvalid);
             %     else
             %         Turb.(ID).wda.(ff{f}) = chi.wda.(ff{f});
             %     end
             % end
         end

         if do_plot
             if ~exist('hfig2', 'var')
                 hfig2 = CreateFigure(is_visible);
                 hfig2.Name = ['Histograms: all final processed estimates'];
             end
             Histograms(Turb.(ID), hfig2, 'pdf', ID, ID);
      
             if do_wda
                 Histograms(Turb.(ID).wda, hfig2, 'pdf', ID, [ID 'W&DA']);
                 hwda = CreateFigure(is_visible);
                 hwda.Name = ['Compare Osborn-Cox vs. Winters-D''Asaro : ' ID];
                 tavg = 3600;
                 ax = plot_estimate(Turb.(ID), ID, tavg);
                 % I want to plot eps calcluated from Kt without
                 % over-complicating plot_estimat; which is set to plot the
                 % eps field.
                 wda_temp = Turb.(ID).wda;
                 wda_temp.eps = wda_temp.eps_Kt;
                 ax = plot_estimate(wda_temp, [ID 'wda'], tavg);
                 legend(ax(3), '\epsilon_\chi', '\epsilon = N^2/\Gamma wda.Kt')
                 % symmetric-log axes for dT/dz and Jq
                 set(ax(5), 'ylim', [-1 1]*2000)
                 symlog(ax(2), 'y', -3);
                 symlog(ax(5), 'y', 1);
                 print(gcf,[basedir '/pics/compare-wda-oc-' ID '.png'],'-dpng','-r200','-painters');
                 if save_fig, savefig(gcf, [basedir '/pics/compare-wda-oc-' ID '.fig']); end
             end

         end
         
         % include statistics (means and medians for each quantity)
         Turb.(ID) = calc_statistics(Turb.(ID));
      end
   end

   if do_plot
       if exist('hfig2', 'var')
          set(0, 'currentfigure', hfig2);
          subplot(221); title(['Final ' num2str(CP.avgwindow/60) ' min mean']);
          subplot(222); title(['Final ' num2str(CP.avgwindow/60) ' min mean']);

          print(gcf,[basedir '/pics/histograms-final.png'],'-dpng','-r200','-painters')
          if save_fig,set(gcf, 'visible', 'on'); savefig(gcf,[basedir '/pics/histograms-final.fig']); end
       end

       if exist('hfraw', 'var')
          set(0, 'currentfigure', hfraw);
          subplot(221); title(['raw 1s estimates']);
          subplot(222); title(['raw 1s estimates']);
          print(gcf,[basedir '/pics/histograms-raw.png'],'-dpng','-r200','-painters')
       end

       if exist('hfstrat', 'var')
          print(hfstrat,[basedir '/pics/histograms-stratification.png'],'-dpng','-r200','-painters')
       end

       if exist('hspecfig', 'var')
           print(hspecfig, [basedir '/pics/histograms-noise-floor.png'],'-dpng','-r200','-painters')
       end
   end

   Turb.hash = githash(['driver' filesep 'combine_turbulence.m']);
   Turb.do_mask = do_mask;
   Turb.parameters = CP;

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
end

%_____________________comparison plot______________________
if do_plot    
    % moved this plotting function to a separate function so it can be run
    % independently if needed.
    [fig] = plot_Turb(basedir, CP.pflag, is_visible);
      print(fig,[basedir '/pics/Compare_Turb.png'],'-dpng','-r200')
      if save_fig
         disp('saving Compare_Turb.fig')
         set(gcf, 'visible', 'on');
         savefig(fig,[basedir '/pics/Compare_Turb.fig'])
      end
end
end

function [CP] = process_estimate_ID(CP, ID)

    if CP.isChipod
        if ID(end) == '1' ...
                    | (length(ID) > 3 & strcmpi(ID(end-3:end), '1_ic')) % sensor T1
            CP.sensor = 1;
        end
        if ID(end) == '2' ...
                    | (length(ID) > 3 & strcmpi(ID(end-3:end), '2_ic')) % sensor T2
            CP.sensor = 2;
        end
    else
        CP.sensor = 1; % gusT has only one sensor
    end

    if ~isempty(strfind(ID, 'p'))
        CP.isPitotEstimate = 1;
    else
        CP.isPitotEstimate = 0;
    end

end

function [chi] = add_nans(CP, chi)
   %___________________NaN out after sensor death______________________
   if CP.sensor == 1 % sensor T1 or gusT T
       death = find(chi.time > CP.T1death, 1, 'first');
       if ~isempty(death)
           disp(['NaNing out sensor T1 after it died on ' datestr(CP.T1death)])
           chi.chi(death:end) = NaN;
           chi.eps(death:end) = NaN;
           chi.T(death:end) = NaN;
           if isfield(chi, 'spec_area')
               chi.spec_area(death:end) = NaN;
           end
       end
       if isfield(chi, 'time_floor')
           deathf = find(chi.time_floor > CP.T1death, 1, 'first');
           if ~isempty(deathf)
               chi.spec_floor(deathf:end) = NaN;
           end
       end
   end

   if CP.isChipod
       if CP.sensor == 2 % sensor T2
           death = find(chi.time > CP.T2death, 1, 'first');
           if ~isempty(death)
               disp(['NaNing out sensor T2 after it died on ' datestr(CP.T2death)])
               chi.chi(death:end) = NaN;
               chi.eps(death:end) = NaN;
               chi.T(death:end) = NaN;
               if isfield(chi, 'spec_area')
                   chi.spec_area(death:end) = NaN;
               end
           end
       end
       if isfield(chi, 'time_floor')
           deathf = find(chi.time_floor > CP.T2death, 1, 'first');
           if ~isempty(deathf)
               chi.spec_floor(deathf:end) = NaN;
           end
       end
   end

   if ~CP.isPitotEstimate
       death = find(chi.time > CP.adcpdeath, 1, 'first');
       if ~isempty(death)
           disp(['NaNing out after mooring velocity died on ' datestr(CP.adcpdeath)]);
           chi.chi(death:end) = NaN;
           chi.eps(death:end) = NaN;
           chi.T(death:end) = NaN;
       end
   end

   %____________________NaN out specific time ranges as necessary____________
   % NaN out spec_area too so that you don't overwrite NaNs with zeros if
   % that condition is satisfied.
   % for temp sensor
   if ~isempty(CP.nantimes{CP.sensor})
       for tt = 1:size(CP.nantimes{CP.sensor}, 1)
           it1 = find_approx(chi.time, CP.nantimes{CP.sensor}(tt, 1), 1);
           it2 = find_approx(chi.time, CP.nantimes{CP.sensor}(tt, 2), 1);
           chi.chi(it1:it2) = NaN;
           chi.eps(it1:it2) = NaN;
           if isfield(chi, 'spec_area')
               chi.spec_area(it1:it2) = NaN;
           end
       end
   end

   if CP.isPitotEstimate & ~isempty(CP.nantimes{3})
       for tt = 1:size(CP.nantimes{3}, 1)
           chi.chi(find_approx(chi.time, CP.nantimes{3}(tt, 1), 1): ...
                   find_approx(chi.time, CP.nantimes{3}(tt, 2), 1)) = NaN;
       end
       if isfield(chi, 'spec_area')
           chi.spec_area(isnan(chi.chi)) = NaN;
       end
       chi.eps(isnan(chi.chi)) = NaN;
   end
end

function [] = compare_min_dTdz_N2_against_SBE_specs(dz, CP, ID)
    % SBE-37 accuracy is 2e-3 C & 3e-3 psu
    sbe_dTdz = 2*2e-3/dz;
    sbe_dSdz = 2*3e-4/dz;
    sbe_N2 = 9.81 * (1.7e-4 * sbe_dTdz + 7.6e-4 * sbe_dSdz);

    if CP.min_dTdz < sbe_dTdz & ID(2) == 'm'
        disp(['WARNING: min_dTdz ' num2str(CP.min_dTdz, '%.1e') ...
              '< minimum resolvable based on ' ...
              'SBE-37 accuracy specifications ' ...
              num2str(sbe_dTdz, '%.1e')]);
    end
    if CP.min_N2 < sbe_N2 & ID(2) == 'm'
        disp(['WARNING: min_N2 ' num2str(CP.min_N2, '%.1e') ...
              ' < minimum resolvable based on ' ...
              'SBE-37 accuracy specifications i.e. ' ...
              num2str(sbe_N2, '%.1e')]);
    end
end

function [chi] = calc_statistics(chi)
    chi.stats.chimean = nanmean(chi.chi);
    chi.stats.epsmean = nanmean(chi.eps);
    chi.stats.Ktmean  = nanmean(chi.Kt);
    chi.stats.Jqmean  = nanmean(chi.Jq);
    chi.stats.chimedian = nanmedian(chi.chi);
    chi.stats.epsmedian = nanmedian(chi.eps);
    chi.stats.Ktmedian  = nanmedian(chi.Kt);
    chi.stats.Jqmedian  = nanmedian(chi.Jq);
    chi.stats.chistd = nanstd(chi.chi);
    chi.stats.epsstd = nanstd(chi.eps);
    chi.stats.Ktstd  = nanstd(chi.Kt);
    chi.stats.Jqstd  = nanstd(chi.Jq);
end

function [Tz] = determine_additional_Tz_mask(CP)

    if CP.additional_mask_dTdz == 'm'
        disp('additional masking using mooring dTdz')
        load([basedir 'input/dTdz_m.mat'])
        Tz = Tz_m;
        clear Tz_m;

    elseif CP.additional_mask_dTdz == 'i'
        disp('additional masking using internal dTdz')
        load([basedir 'input/dTdz_i.mat'])
        Tz = Tz_i;
        clear Tz_i;
    else
        Tz = [];
    end

end

function [spdmask, addspdmask] = determine_speed_masks(basedir, ID, CP, chi)

    % speed mask could change depending on estimate
    mask_spd = ID(1);

    if mask_spd == 'm'
        load([basedir '/input/vel_m.mat']);
        vel = vel_m;
    elseif mask_spd == 'p' ...
            & exist([basedir '/input/vel_p.mat'], 'file')
        load([basedir '/input/vel_p.mat']);
        vel = vel_p;
    else
        % hack to avoid background velocity masking
        % if vel_p.mat does not exist
        vel.time = chi.time;
        vel.spd = ones(size(chi.time));
    end
    spdmask = interp1(vel.time, vel.spd, chi.time);

    if CP.additional_mask_spd ~= ''
        if CP.additional_mask_spd == 'm' & ~exist('vel_m', 'var')
            load([basedir '/input/vel_m.mat']);
            vel = vel_m;
        elseif CP.additional_mask_spd == 'p' & ~exist('vel_p', 'var')
            load([basedir '/input/vel_p.mat']);
            vel = vel_p;
        end
        addspdmask = interp1(vel.time, vel.spd, chi.time);
    else
        addspdmask = [];
    end
end

function [motion] = process_motion_file(motionfile)
    if ~exist(motionfile, 'file')
        error(['proc/motion.mat not found. IGNORING ' ...
               'FLUSHING MASKING. Run do_temp_proc ' ...
               'to create proc/motion.mat']);
    end

    if ~isequal(motion.time, chi.time)
        ff = fieldnames(motion);
        for f=1:length(ff)
            if strcmpi(ff{f}, 'time'), continue; end
            motion.(ff{f}) = interp1(motion.time, motion.(ff{f}), chi.time);
        end
        motion.time = chi.time;
    end
end

function [out] = truncate_time(in, time_range)

    % find desired time range
    iiTrange = find( in.time >= time_range(1) & in.time <= time_range(2) );

    % extract valid time range here so that histograms work
    ff = fieldnames(in);
    for nn = 1:length(ff)
        try
            out.(ff{nn}) = in.(ff{nn})(iiTrange);
        catch ME;
            out.(ff{nn}) = in.(ff{nn});
        end
    end
    
    if isfield(in,'time_floor')
        iifloor = find( in.time_floor >= time_range(1) & in.time_floor <= time_range(2) );
        out.time_floor = in.time_floor(iifloor);
        out.spec_floor = in.spec_floor(iifloor);
    end

end

function [] = make_noise_floor_histograms(ID, chi, hspec1, hspec2, hspec3, ...
                                          spec_floor, spec_floor_mask, do_spec_area)

    clr = choose_color(ID, 'color');

    if do_spec_area
        histogram(hspec1, log10(chi.spec_area), 'EdgeColor', clr, ...
                  'displaystyle', 'stairs', 'DisplayName', ID);
        plot(hspec1, [1 1]* log10(nanmedian(spec_floor) * nanmean(chi.nfft)), hspec1.YLim, ...
             '--', 'color', clr, 'displayname', [ID ' noise floor'])
        legend(hspec1, '-dynamiclegend');
        xlabel(hspec1, 'log_{10} area under Tp spectrum')
        ylabel(hspec1, 'count')
    end

    histogram(hspec2, log10(chi.chi(~spec_floor_mask)), 'EdgeColor', clr, ...
              'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(hspec2, log10(chi.chi(spec_floor_mask)), 'linestyle', '--', ...
              'EdgeColor', clr, 'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(hspec2, log10(chi.chi), 'EdgeColor', clr, 'linewidth', 2, ...
              'normalization', 'count', 'displaystyle', 'stairs', 'displayname', ['\chi_{' ID '}'])
    legend(hspec2, '-dynamiclegend')
    xlabel(hspec2, 'log_{10} \chi')
    ylabel(hspec2, 'count')
    title(hspec2, 'thin dashed, solid = (below, above) noise floor | thick = all values')

    histogram(hspec3, log10(chi.eps(~spec_floor_mask)), 'EdgeColor', clr, ...
              'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(hspec3, log10(chi.eps(spec_floor_mask)), 'linestyle', '--', ...
              'EdgeColor', clr, 'normalization', 'count', 'displaystyle', 'stairs', ...
              'HandleVisibility', 'off')
    histogram(hspec3, log10(chi.eps), 'EdgeColor', clr, 'linewidth', 2, ...
              'normalization', 'count', 'displaystyle', 'stairs', 'displayname', ['\epsilon_{' ID '}'])
    legend(hspec3, '-dynamiclegend')
    xlabel(hspec3, 'log_{10} \epsilon')
    ylabel(hspec3, 'count')

end

% Examples of using TestMask and plot_estimate to check masking
% TestMask(chi, abs(chi.dTdz), '<', [1e-4, 3e-4, 1e-3], 'Tz');
% t0 = datenum(2014, 01, 01);
% t1 = datenum(2014, 03, 01);
% plot_estimate(chi, 'raw', 0, [], t0, t1);
% chi1 = ApplyMask(chi, abs(chi.dTdz), '<', 3e-4, 'T_z');
% plot_estimate(chi1, 'T_z > 3e-4', 0, [], t0, t1);
% load ../proc/temp.mat
% load ../proc/Turb.mat
