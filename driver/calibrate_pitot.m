%     this is the driver for Pitot calibration
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 11:30:20 PST 2017

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_raw_data = 0;     % do the averaging of the raw-data (1) or skip (0) if done before
   do_v0_self  = 0;     % detremine V0 based on a min of the averaged signal (self contained)
   do_v0_adcp  = 0;     % detremin V0 based on a fit against reference velocity (adcp) data
   do_plot     = 0;     % generate some figures in ../pics/ to compare the different velocity estimates

   % This is the time range where the pitot sensor is returning
   % good data
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0);
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0);

   % calibrate in time range different from valid data time range?
   % if so set limits here just as for time_range.
   % by default, both time ranges are equal.
   cal_time_range = time_range;

   % which temperature sensor to use T1 (1) or if T1 is broken T2 (2) ;  
   % for gusTs (0)
   use_T = 1;  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%_____________________include path of processing files______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   
%_____________________get header files______________________
   % general headers
      hfid = [basedir filesep 'calib' filesep 'header.mat'] ;
      if exist(hfid, 'file')
         load(hfid)
      else % no header found
         error(['I could not load ' hfid ' !!!!Run create_header.m first !!!!']);
      end

   % Pitot header
      hfid = [basedir filesep 'calib' filesep 'header_p.mat'] ;
      if exist(hfid, 'file')
         load(hfid)
      else % no header found
         error(['I could not load ' hfid ' !!!!Run create_header.m first !!!!']);
      end

      % Make a backup of old Pitot header file
      if isfield(W, 'V0' )
         disp('The header_p file does contain already V0');
         disp('The header_p file is saved as a back up in');
         disp('./calib/header_p.mat.backup');
         disp('... And we recalculate V0');
         save([hfid '.backup'], 'W');
      end

      if do_raw_data
          do_raw_pitot;
      end
%_____________________load averaged raw data and do basic calibration______________________
   fid = ['../proc/Praw.mat'];
   if exist(fid, 'file');
      load(fid);
   else
      disp(['The raw data are not processed yet you need to set']);
      disp('do_raw_data = 1; or run do_raw_pitot.m')
   end
   %--------------------base calibation----------------------
   P.time = Praw.time;

   % temperature
   if (use_T==0) %gusTs
      P.T   =  (Praw.T.^2+Praw.vT)*head.coef.T(3)+ Praw.T*head.coef.T(2) + head.coef.T(1);
   else  % chipods
      P.T1   =  (Praw.T1.^2+Praw.vT1)*head.coef.T1(3)+ Praw.T1*head.coef.T1(2) + head.coef.T1(1);
      P.T2   =  (Praw.T2.^2+Praw.vT2)*head.coef.T2(3)+ Praw.T2*head.coef.T2(2) + head.coef.T2(1);
      if use_T==1
         P.T    = P.T1; 
      else % in case T1 is broken
         P.T    = P.T2;
      end
   end

   % pressure
   P.P    =  Praw.P*head.coef.P(2) + head.coef.P(1);

   % compass
      if isfield(head.coef, 'CMP')
         P.cmp  = Praw.cmp + head.coef.CMP(1);
      else
         P.cmp  = Praw.cmp;
         disp(['CMP' ' does not exit in header']);
      end

   %---------------------pre calibration for Pitot----------------------
      %% find all idexes in the desired time interval;
      iiPcal = find( P.time>=cal_time_range(1) & P.time<=cal_time_range(2) );
      iiP = find( P.time>=time_range(1) & P.time<=time_range(2) );

   % set the average temperature as reference value for the Pitot calibration
   W.T0   =  nanmean(P.T(iiPcal));
   W.P0   =  nanmean(P.P(iiPcal));

   % calibrate the Pitot voltage for temperature (pressure ? Tilt ?)
   P.W   =   Praw.W - (P.T-W.T0)*W.T(2);



%_____________________detremine V0 based on min method (self contained)______________________
Porg = P;
if do_v0_self

   % calculate V0 as the median of the smallest 5 % of the averaged values
      w_sort = sort(P.W(iiPcal));
      W.V0 = median(w_sort(1:round(length(w_sort)/20)));

   % calibrate voltage into speeds
   [P.spd, ~, ~] = pitot_calibrate(P.W, P.T, 0, W.V0, W.T0, 0, 1/W.Pd(2), 0, 0);

   % add directional information from the compass
   P.U = pitot_add_direction(P.time, P.spd, P.time, P.cmp);


   % output
   disp(['based on the internal method V0 is calculated to be']);
   W
   
   if do_plot
       CreateFigure;
         a=1;
         ax(a) = subplot(3,1,a);
            plot(ax(a), P.time, P.T, 'Linewidth', 1);
            hold all;
            plot(ax(a), P.time([1 end]), [1 1]*W.T0, 'Linewidth', 1);
            ylabel(ax(a), 'T [deg C]');
            datetick(ax(a), 'keeplimits');
            legend(ax(a),  'T signal', 'T_0');
         a=2;
         ax(a) = subplot(3,1,a);
            plot(ax(a), P.time, P.P, 'Linewidth', 1);
            hold all;
            plot(ax(a), P.time([1 end]), [1 1]*W.P0, 'Linewidth', 1);
            ylabel(ax(a), 'Pres [psu]');
            legend(ax(a),  'P signal', 'P_0');
            datetick(ax(a), 'keeplimits');
         a=3;
         ax(a) = subplot(3,1,a);
            plot(ax(a), P.time, P.W, 'Linewidth', 1);
            hold all;
            plot(ax(a), P.time([1 end]), [1 1]*W.V0, 'Linewidth', 1);
            ylabel(ax(a), '[Volt]');
            legend(ax(a),  'Pitot signal', 'V_0');
            datetick(ax(a), 'keeplimits');

            linkaxes(ax, 'x');
            xlim(ax(1), time_range)
            
   end

   % cut data matrix
   ff = fields(P);
   for fi = 1:length(ff)
      P.(ff{fi})   = P.(ff{fi})(iiP);
   end

   % save header and calibrated data
   save('../calib/header_p_self.mat', 'W');
   save('../calib/header_p.mat', 'W');

   save('../proc/P_self.mat', 'P');
   P = Porg;
end

%_____________________detremine V0 based on a fit against ADCP data______________________
if do_v0_adcp
  
   % load vel_m.mat as reference velocity data
      load ../input/vel_m.mat;
   

   %% find all idexes in the desired time interval;
      % adcp
      iiA = find( vel_m.time>=P.time(iiPcal(1)) &  vel_m.time<=P.time(iiPcal(end)) );

   % determine V0
       [W.V0] = fit_pitot_v0( vel_m.time, vel_m.spd, P.time(iiPcal), P.W(iiPcal), 1/W.Pd(2), do_plot);
       if do_plot, print(gcf, '../pics/pitot-adcp-fit-voltages.png', '-dpng', '-r200', '-painters'); end
   % calibrate voltage into speeds
   [P.spd, ~, ~] = pitot_calibrate(P.W, P.T, 0, W.V0, W.T0, 0, 1/W.Pd(2), 0, 0);

   % add directional information from the compass
   P.U = pitot_add_direction(P.time, P.spd, P.time, P.cmp);

   % cut data matrix
   ff = fields(P);
   for fi = 1:length(ff)
      P.(ff{fi})   = P.(ff{fi})(iiP);
   end

   % save header and calibrated data
   save('../calib/header_p_fit.mat', 'W');
   save('../calib/header_p.mat', 'W');

   save('../proc/P_fit.mat', 'P');

   % output
   disp(['based on the fitting method V0 is calculated to be']);
   W


   % calculate direction off-set
   D_off = (angle(nanmean(P.U)) - angle(nanmean(vel_m.U(iiA))))*180/pi;
   disp('The direction off set between ADCP and Chipod is');
   disp([num2str( D_off ) ' deg']);

   if do_plot
      % generate a comparison plot
      a_L    = 'ADCP';
      p_L    = 'Pitot V0_{fit}';
      [fig] =  compare_velocity_timeseries(vel_m.time, vel_m.U, a_L, P.time, P.U, p_L);
      print(gcf,'../pics/Pitot_vs_ADCP_V0_fit.png','-dpng','-r200','-painters')
   end

end


%_____________________compare differnt methods______________________
if (do_v0_adcp & do_v0_self & do_plot)

   Ps = load('../proc/P_self.mat');
   Pf = load('../proc/P_fit.mat');

      a_L    = 'V0_{fit}';
      p_L    = 'V0_{self}';
      [fig] =  compare_velocity_timeseries(Pf.P.time, Pf.P.U, a_L, Ps.P.time, Ps.P.U, p_L);
      print(gcf,'../pics/V0_fit_vs_self.png','-dpng','-r200','-painters')
end

%___________________generating Pitot velocity input_________________
fidf = '../proc/P_fit.mat';
fids = '../proc/P_self.mat';

if exist(fidf, 'file');
    load(fidf);
    vel_p.text = 'vel_p.mat is generated based on the ADCP fitted Pitot signal';
    disp(vel_p.text);
elseif exist(fids, 'file');
    load(fids);
    vel_p.text = 'vel_p.mat is generated in the self contained way';
    disp(vel_p.text);
else
    error(['Pitot velocities have not been calibrated yet!. Run with do_v0_adcp = 1 ' ...
           'or do_v0_self = 1 first!']);
end

vel_p.time  = P.time;
vel_p.spd   = P.spd;
vel_p.U     = P.U;
vel_p.u     = real(P.U);
vel_p.v     = imag(P.U);

save('../input/vel_p.mat', 'vel_p');
disp('vel_p.mat created!')