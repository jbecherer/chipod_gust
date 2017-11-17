function [data] = chi_calibrate_chipod(rfid, head)
%%
%     
%        This function reads CHIPOD raw data and calibrates them acording
%        to the coefficients in head.
%
%        Input
%           rfid   :  path the the specific raw-file
%           head   :  corresponding header
%
%        Output
%           data   : data structure containing calibrated data
%
%   created by: 
%        Johannes Becherer
%        Fri Sep  2 15:53:26 PDT 2016


%_____________________read in raw data______________________

[rdat, headtest]  = raw_load_chipod(rfid);

%____________ check if raw_file head corresponds to given head___________

% if ~isequal(head.coef, headtest.coef)
%    disp(' ');
%    disp('!!! WARNING !!!!!');
%    disp('the raw file header is different from the given header');
%    disp('the given header is used and the raw-file header ignored');
%    disp(' ');
% end

%_____________________Syncronize time______________________
   % data on Tp time stemp (datenum, T1P, T2P)
   %         T1 time stemp (T1, T2, AX, AY, AZ, W, WP. P)
   %         CMP time stemp (CMP)

   % time vector
      chi.time_tp  = rdat.datenum;
         % remove bad time glitches
         Ntp   = length(chi.time_tp); 
         Itp   = 1:Ntp;
         Ibad  = ( chi.time_tp<(nanmedian(chi.time_tp)-1) |...
                        chi.time_tp>(nanmedian(chi.time_tp)+1) | isnan(chi.time_tp) );    % usually time glitches are more than one day off
         if ~isempty(Ibad)
            chi.time_tp(Ibad) =  interp1( Itp(~Ibad), chi.time_tp(~Ibad), Itp(Ibad), 'linear','extrap');
         end
         % glitch data are also bad in the Tp signal 
         rdat.T1P(Ibad) = nan; 
         rdat.T2P(Ibad) = nan; 
         % other signals are also glitched with 0 values 

      % syncronize time by making Tp time step must be divdable by 20 to get syncronized times 
      chi.time_tp  = chi.time_tp(1:(floor(Ntp/20)*20));
         Ntp   = length(chi.time_tp); 
      chi.time     = chi.time_tp(1:2:end);
         Nt = length(chi.time);
      chi.time_cmp = chi.time(1:10:end);
         Ntc = length(chi.time_cmp);
%_____________________calibrate data______________________


   % tempertuare
      chi.T1=calibrate_polynomial(rdat.T1,head.coef.T1);
      chi.T2=calibrate_polynomial(rdat.T2,head.coef.T2);
         % remove glitches
         chi.T1(chi.T1==0) = nan;
         chi.T2(chi.T2==0) = nan;
         % time sync
         chi.T1 = chi.T1(1:Nt);
         chi.T2 = chi.T2(1:Nt);
   % pressure
      chi.P=calibrate_polynomial(rdat.P,head.coef.P);
         % remove glitches
         chi.P(chi.P==0) = nan;
         % time sync
         chi.P = chi.P(1:Nt);

          % integrate pressure sensor
          [chi.p_dis_z, chi.p_vel_z] = integrate_pres(chi, head);

      
   % accelerometer
         g=9.81;
         chi.AX=g.*calibrate_polynomial(rdat.AX,head.coef.AX);
         chi.AY=g.*calibrate_polynomial(rdat.AY,head.coef.AY);
         chi.AZ=g.*calibrate_polynomial(rdat.AZ,head.coef.AZ);
            % remove glitched data
            chi.AX(chi.AX==0) = nan;
            chi.AY(chi.AY==0) = nan;
            chi.AZ(chi.AZ==0) = nan;
             chi.AX=fillgap(chi.AX);
             chi.AY=fillgap(chi.AY);
             chi.AZ=fillgap(chi.AZ);
            % time sync
            chi.AX = chi.AX(1:Nt);
            chi.AY = chi.AY(1:Nt);
            chi.AZ = chi.AZ(1:Nt);
          [dis,vel]=integrate_acc(chi,head);

          chi.Acc = sqrt((chi.AX-nanmean(chi.AX)).^2 + (chi.AY-nanmean(chi.AY)).^2 + (chi.AZ-nanmean(chi.AZ)).^2 );

          chi.a_dis_x = dis.x;
          chi.a_dis_y = dis.y;
          chi.a_dis_z = dis.z;
          chi.a_vel_x = vel.x;
          chi.a_vel_y = vel.y;
          chi.a_vel_z = vel.z;
           
         chi.AXtilt=calibrate_tilt(rdat.AX,head.coef.AX);
         chi.AYtilt=calibrate_tilt(rdat.AY,head.coef.AY);
         chi.AZtilt=calibrate_tilt(rdat.AZ,head.coef.AZ);

          % compare AZ with dP/dt
          % figure;
          % subplot(121)
          % histogram2(chi.a_vel_z, chi.p_vel_z, 'displaystyle', 'tile', ...
          %            'normalization', 'pdf')
          % axis square;
          % xlabel('accel vel_z'); ylabel('pres vel_z')
          % line45
          % title(['PDF | r = ' ...
          %        num2str(min(min(corrcoef(chi.a_vel_z, chi.p_vel_z))))])
          % subplot(122)
          % histogram2(chi.a_dis_z, chi.p_dis_z, 'displaystyle', ...
          %            'tile', 'normalization', 'pdf')
          % xlabel('accel disp_z'); ylabel('pres disp_z')
          % title(['PDF | r = ' ...
          %        num2str(min(min(corrcoef(chi.a_dis_z, chi.p_dis_z))))])
          % axis square;
          % line45

         % compass
         chi.cmp = rdat.CMP/10+head.coef.CMP(1);
         % In case time_cmp is to long (sometimes its one index to long)
         if length(chi.cmp)< Ntc
            chi.time_cmp = chi.time_cmp(1:length(chi.cmp));
         elseif length(chi.cmp)> Ntc
            chi.cmp = chi.cmp(1:Ntc);
         end
      

   % DTdt
      rdat.T1P = rdat.T1P - nanmean(rdat.T1P);
      rdat.T2P = rdat.T2P - nanmean(rdat.T2P);
         % sync time
         rdat.T1P = rdat.T1P(1:Ntp); 
         rdat.T2P = rdat.T2P(1:Ntp); 
         rdat.T1 = rdat.T1(1:Nt); 
         rdat.T2 = rdat.T2(1:Nt); 
      chi.T1Pt = calibrate_tp( rdat.T1P, head.coef.T1P, rdat.T1, head.coef.T1, 100*ones(size(rdat.T1)) );
      chi.T2Pt = calibrate_tp( rdat.T2P, head.coef.T2P, rdat.T2, head.coef.T2, 100*ones(size(rdat.T1)) );

   % pitot_voltage
      % find pitot data W or WP
       if isfield(rdat, 'W')
         dV1 = abs(nanmean(rdat.W)-2.02);
         dV2 = abs(nanmean(rdat.WP)-2.02);
         if dV1>dV2
            chi.W  = rdat.W;
         else
            chi.W  = rdat.WP;
         end
       else  
         dV1 = abs(nanmean(rdat.W2)-2.02);
         dV2 = abs(nanmean(rdat.W3)-2.02);
         if dV1>dV2
            chi.W  = rdat.W2;
         else
            chi.W  = rdat.W3;
         end
       end

       % sync time
       chi.W  = chi.W(1:Nt);

% change unit of pressure
chi.depth = (chi.P-14.7)/1.47;
%---------------------return data----------------------
data = chi;


