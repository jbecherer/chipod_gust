%     this is the driver for Pitot calibration
%
%   created by: 
%        Johannes Becherer
%        Tue Feb 14 11:30:20 PST 2017

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_raw_data = 0;     % do the averagring of the raw-data (1) or skip (0) if done before 
   do_v0_self  = 0;     % detremine V0 based on a min of the averaged signal (self contained)
   do_v0_adcp  = 0;     % detremin V0 based on a fit against reference velocity (adcp) data
   do_plot     = 0;     % generate some figures in ../pics/ to compare the different velocity estimates

   % if you want to restrict the time range that should be analyzed use the following
   time_range(1)  = datenum(2000, 1, 1, 0, 0, 0); 
   time_range(2)  = datenum(2030, 1, 1, 0, 0, 0); 


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
      
%_____________________average all Pitot raw data to 10 min bins______________________

if do_raw_data

      disp('averaging Pitot raw  data in ./proc/Praw.mat ')
      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
           try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);

               if use_T == 0 % gusT
                  [~] = pitot_avg_raw_data(basedir, fids{f}, 0);
               else        %chipod
                  [~] = pitot_avg_raw_data(basedir, fids{f});
               end
           catch
               disp(['!!!!!! ' fids{f} ' crashed while processing vel_p structure !!!!!!' ]);
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               if use_T == 0 % gusT
                  [Praw] = pitot_avg_raw_data(basedir, fids{f}, 0);
               else        %chipod
                  [Praw] = pitot_avg_raw_data(basedir, fids{f});
               end
         end
      end
   %_____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'Praw', 0);
end

%_____________________load averaged raw data and do basic calibration______________________
   fid = ['../proc/Praw.mat'];
   if exist(fid, 'file');
      load(fid);
   else
      disp(['The raw data are not processed yet you need to set']);
      disp('do_raw_data = 1;')
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
      iiP = find( P.time>=time_range(1) & P.time<=time_range(2) );

   % set the average temperature as reference value for the Pitot calibration
   W.T0   =  nanmean(P.T(iiP));
   W.P0   =  nanmean(P.P(iiP));

   % calibrate the Pitot voltage for temperature (pressure ? Tilt ?)
   P.W   =   Praw.W - (P.T-W.T0)*W.T(2);



%_____________________detremine V0 based on min method (self contained)______________________
Porg = P;
if do_v0_self

   % calculate V0 as the median of the smallest 5 % of the averaged values
      w_sort = sort(P.W(iiP));
      W.V0 = median(w_sort(1:round(length(w_sort)/20)));

   % calibrate voltage into speeds
   % temperature calibration done earlier so set that to 0
   [P.spd, ~, ~] = pitot_calibrate(P.W, P.T, 0, W.V0, W.T0, 0, 1/W.Pd(2), 0, 0);

   disp(['Time instants with calibrated Pd < 0 = ' ...
        num2str(sum(P.spd(iiP) == 0)/length(P.spd(iiP))*100) '%'])
   % add directional information from the compass
   P.U = pitot_add_direction(P.time, P.spd, P.time, P.cmp);


   % output
   disp(['based on the internal method V0 is calculated to be']);
   W
   
   if do_plot
      figure
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
      iiA = find( vel_m.time>=P.time(iiP(1)) &  vel_m.time<=P.time(iiP(end)) );


   % determine V0
       [W.V0] = fit_pitot_v0( vel_m.time, vel_m.spd, P.time(iiP), P.W(iiP), 1/W.Pd(2), 0);
   
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

   if abs(D_off)>2 % Ask if the chipod header should be changed
         choice = questdlg(['I found a direction off set between chipod and ADCP of ' num2str( D_off ) ' deg.' ...
                           ' Do you want me to change the compass off set in the chipod header?'], ...
                           'Compass off set', ...
                            'Yes','No','No');
         switch choice 
            case 'Yes'
               head.coef.CMP(1) = head.coef.CMP(1) + D_off;
               save('../calib/header.mat', 'head');
               disp('Chipod header changed!')
            case 'No'
               disp('Chipod header NOT changed!')
         end
               
   end

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
