function [] = determine_v0( basedir, do_v0_self, do_v0_adcp,  do_plot, do_vel_p, time_range, use_T, use_press, vis, DcalWindow, DcalIncrement )
%% [] = determine_v0( basedir, do_v0_self, do_v0_adcp, do_plot, do_vel_p, time_range )
%     
%     This function is meant to determine V0 the pitot voltage off set based on different methods
%
%
%  INPUT
%     basedir     :   base directory
%     do_v0_self  :   detremine V0 based on a min of the averaged signal (self contained)
%     do_v0_adcp  :   detremin V0 based on a fit against reference velocity (adcp) data 
%     do_plot     :   generate some figures in ../pics/ to compare the different velocity estimates
%     do_vel_p    :   which calibration should be used for vel_p (0 none (default), 1: adcp, 2: self)
%     use_T       :   which T sensor should be used (default 1)
%     use_press   :   shall the pressure calibration be used (default 1)
%     vis         :  shall the figures be vissble (default 'on', 'off')
%
%
%   created by: 
%        Johannes Becherer
%        Tue Nov 28 14:26:43 PST 2017


if nargin < 7
    use_T = 1;
end
if nargin < 8
    use_press = 1;
end
if nargin < 9
    vis = 'on';
end
if nargin < 10
     DcalWindow =  1000;    % 1000 days window for V0 time (effectivly only a single value for the entire record)
end
if nargin < 11
     DcalIncrement = 1000;  % 1000 day increment
end

fidf    = [basedir '/proc/P_fit.mat'];
fids    = [basedir '/proc/P_self.mat'];
fidvelm = [basedir '/input/vel_m.mat'];
if ~exist(fidvelm)
   do_v0_adcp = 0;
   disp('vel_m.mat does not exist there for do_v0_adcp = 0');
end

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
      elseif isfield(head.coef, 'W') %maybe header for pitot hidden in normal header 
         W = head.coef.W;
         warning('Could not find header_p.mat. Used head.coef.W instead!')
      else  % no header found
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

      % remove pressure calibration
      if ~use_press & W.Ps(2) ~=0
         disp('!!! the pressure calibration is switched off by setting W.Ps(2) = 0')
         W.Ps(2)  =  0;
      end

%_____________________load averaged raw data and do basic calibration______________________
   fid = [basedir '/proc/Praw.mat'];
   if exist(fid, 'file');
      load(fid);
   else
      disp(['The raw data are not processed yet you need to set']);
      disp('do_raw_data = 1; or run do_raw_pitot.m')
   end





   %--------------------base calibation----------------------
   if time_range(1)< Praw.time(1)
      time_range(1) =Praw.time(1);  
   end
   if time_range(2)>  Praw.time(end)
      time_range(2) =Praw.time(end);  
   end

   ii_time_range  =  find( Praw.time>=time_range(1) & Praw.time<=time_range(2) );
   P.time         = Praw.time(ii_time_range);

   % temperature
   if (isfield(Praw, 'T')) %gusTs
      P.T   =  (Praw.T(ii_time_range).^2+Praw.vT(ii_time_range))*head.coef.T(3)+ Praw.T(ii_time_range)*head.coef.T(2) + head.coef.T(1);
   else  % chipods
      P.T1   =  (Praw.T1(ii_time_range).^2+Praw.vT1(ii_time_range))*head.coef.T1(3)+ Praw.T1(ii_time_range)*head.coef.T1(2) + head.coef.T1(1);
      P.T2   =  (Praw.T2(ii_time_range).^2+Praw.vT2(ii_time_range))*head.coef.T2(3)+ Praw.T2(ii_time_range)*head.coef.T2(2) + head.coef.T2(1);
      if use_T==1
         P.T    = P.T1; 
      else % in case T1 is broken
         P.T    = P.T2;
      end
   end

   %% recalculate slope for reduced temperature range
   if isfield(W, 'data')
      ii_T_range  =  find( W.data.T.degC>=min(P.T) & W.data.T.degC<=max(P.T) );
      T_coeffs_new   =  polyfit( W.data.T.degC, W.data.T.V, 1  );
      disp([' Calculated new temperature slope : ' num2str(T_coeffs_new(1))]);
      disp(['                        old slope : ' num2str(W.T(2))]);
      W.T(2) = T_coeffs_new(1);
   end

   % pressure
   P.P    =  Praw.P(ii_time_range)*head.coef.P(2) + head.coef.P(1);

   % compass
      if isfield(head.coef, 'CMP')
         P.cmp  = Praw.cmp(ii_time_range) + head.coef.CMP(1);
      else
         P.cmp  = Praw.cmp(ii_time_range);
         disp(['CMP' ' does not exit in header']);
      end

   %---------------------pre calibration for Pitot----------------------
      %% find all idexes in the desired time interval;
      iiPcal = find( P.time>=time_range(1) & P.time<=time_range(2) );
      iiP = find( P.time>=time_range(1) & P.time<=time_range(2) );

   % set the average temperature as reference value for the Pitot calibration
   W.T0   =  nanmean(P.T(iiPcal));
   W.P0   =  nanmean(P.P(iiPcal));

   % calibrate the Pitot voltage for temperature (pressure ? Tilt ?)
   P.W   =   Praw.W(ii_time_range) - (P.T-W.T0)*W.T(2);


%_____________________self methode______________________
Porg = P;

P.spd = nan(size(P.W));
if do_v0_self
      W.V0     =  nan; 
      W.time   =  nan; 
      
      cnt   =  1;
      time_low =  time_range(1);
      time_up  =  time_low+ DcalWindow;
         % for an interval larger than the record
         if time_up > time_range(2)
            time_up = time_range(2)
         end
      while time_up <= time_range(2)
         iiWcal   =  find( P.time>=time_low & P.time<=time_up );
         if ~isempty(iiWcal)
         W.time(cnt)   =  nanmean(P.time(iiWcal));
         W.V0(cnt)  =  v0_self(P.W(iiWcal));

         % calibrate voltage into speeds
         % temperature calibration done earlier so set that to 0
         W1 = W;
            W1.P0 = 0; % switch off temp and press calibration
            W1.T = [0 0 0 0 0];
            W1.Ps = [0 0 0 0 0];
            W1.V0 =  W.V0(cnt);
         [P.spd(iiWcal), ~, ~] = pitot_calibrate(P.W(iiWcal), P.T(iiWcal), 0, W1);

         cnt = cnt+1;
         end
         time_low =  time_low + DcalIncrement;
         time_up  =  time_low + DcalWindow;
      end


   disp(['Time instants with calibrated Pd < 0 = ' ...
        num2str(sum(P.spd(iiP) == 0)/length(P.spd(iiP))*100) '%'])
   % add directional information from the compass
   P.U = pitot_add_direction(P.time, P.spd, P.time, P.cmp);


   % output
   disp(['based on the internal method V0 is calculated to be']);
   W
   
   if do_plot
       CreateFigure(vis);
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
            plot(ax(a), W.time, W.V0, 'Linewidth', 1);
            plot(ax(a), W.time, W.V0,'+', 'Linewidth', 1);
            if length(W.V0) == 1
               plot(ax(a), time_range, [1 1]*W.V0, 'Linewidth', 1);
            end

            ylabel(ax(a), '[Volt]');
            legend(ax(a),  'Pitot signal', 'V_0');
            datetick(ax(a), 'keeplimits');

            linkaxes(ax, 'x');
            xlim(ax(1), time_range)

         print(gcf,[basedir '/pics/pitot_self_diagnostic.png'],'-dpng','-r100','-painters')
         
            
   end

   % cut data matrix
   ff = fields(P);
   for fi = 1:length(ff)
      P.(ff{fi})   = P.(ff{fi})(iiP);
   end

   % save header and calibrated data
   save([basedir '/calib/header_p_self.mat'], 'W');
   save([basedir '/calib/header_p.mat'], 'W');

   save([basedir '/proc/P_self.mat'], 'P');
   P = Porg;
end

%_____________________detremine V0 based on a fit against ADCP data______________________
if do_v0_adcp
  
   % load vel_m.mat as reference velocity data
      if exist(fidvelm)
         load(fidvelm);
      else
         error([fidvelm ' does not exist'])
      end
   

   %% find all idexes in the desired time interval;
      % adcp
      iiA = find( vel_m.time>=P.time(iiPcal(1)) &  vel_m.time<=P.time(iiPcal(end)) );

   % determine V0
       [W.V0] = fit_pitot_v0( vel_m.time, vel_m.spd, P.time(iiPcal), P.W(iiPcal), 1/W.Pd(2), do_plot, vis);
       if do_plot, print(gcf, [basedir '/pics/pitot-adcp-fit-voltages.png'], '-dpng', '-r200'); end
   % calibrate voltage into speeds
   W1 = W;
      W1.P0 = 0; % switch off temp and press calibration
      W1.T = [0 0 0 0 0];
      W1.Ps = [0 0 0 0 0];
   [P.spd, ~, ~] = pitot_calibrate(P.W, P.T, 0, W1);

   % add directional information from the compass
   P.U = pitot_add_direction(P.time, P.spd, P.time, P.cmp);

   % cut data matrix
   ff = fields(P);
   for fi = 1:length(ff)
      P.(ff{fi})   = P.(ff{fi})(iiP);
   end

   % save header and calibrated data
   save([basedir '/calib/header_p_fit.mat'], 'W');
   if ~exist([basedir '/calib/header_p_self.mat']) % use only ADCP-fit header if self is not available
      save([basedir '/calib/header_p.mat'], 'W');
   end

   save([basedir '/proc/P_fit.mat'], 'P');

   % output
   disp(['based on the fitting method V0 is calculated to be']);
   W


   % calculate direction off-set
   D_off = (angle(nanmean(P.U)) - angle(nanmean(vel_m.U(iiA))))*180/pi;
   disp('The direction off set between ADCP and Chipod is');
   disp([num2str( D_off ) ' deg']);


end


%_____________________compare differnt methods______________________
if do_plot

   if (exist(fidf, 'file') & exist(fids, 'file') )

      load('../input/vel_m.mat');
      Ps = load(fids);
      Pf = load(fidf);

         a_L    = 'V0_{fit}';
         p_L    = 'V0_{self}';
         [fig] =  compare_velocity_timeseries(Pf.P.time, Pf.P.U, a_L, Ps.P.time, Ps.P.U, p_L, vis);
         print(gcf,[basedir '/pics/pitot_V0_fit_vs_self.png'],'-dpng','-r100')
         %savefig(fig,[basedir '/pics/pitot_V0_fit_vs_self.fig'])
      a_L    = 'ADCP';
      p_L    = 'Pitot V0_{self}';
      [fig] =  compare_velocity_timeseries(vel_m.time, vel_m.U, a_L, Ps.P.time, Ps.P.U, p_L, vis);
      print(fig,[basedir '/pics/pitot_self_vs_ADCP_V0_fit.png'],'-dpng','-r100')
   end
   if (exist(fidf, 'file') & exist(fidvelm, 'file') )
      load(fidvelm);
      % generate a comparison plot
      Pf = load([basedir '/proc/P_fit.mat']);
      a_L    = 'ADCP';
      p_L    = 'Pitot V0_{fit}';
      [fig] =  compare_velocity_timeseries(vel_m.time, vel_m.U, a_L, Pf.P.time, Pf.P.U, p_L, vis);
      print(fig,[basedir '/pics/pitot_fit_vs_ADCP_V0_fit.png'],'-dpng','-r100')
      %savefig(fig,[basedir '/pics/pitot_vs_ADCP_V0_fit.fig'])
   end
end

%___________________generating vel_p________________
if do_vel_p > 0
   if ( do_vel_p == 2 & exist(fids, 'file'));
       load(fids);
       vel_p.text = 'vel_p.mat is generated in the self contained way';
       disp(vel_p.text);
   elseif ( do_vel_p == 1 & exist(fidf, 'file'));
       load(fidf);
       vel_p.text = 'vel_p.mat is generated based on the ADCP fitted Pitot signal';
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

   save([basedir '/input/vel_p.mat'], 'vel_p');
   disp('vel_p.mat created!')
end


end % end of main function


function [V0]  =  v0_self(W)
   % calculate V0 as the median of the smallest 5 % of the averaged values
      w_sort = sort(W);
      V0 = median(w_sort(1:round(length(w_sort)/20)));
end

