function [TL] =   whoAmI_timeLimits(basedir)
%%    [TL] =   whoAmI_timeLimits(basedir)
%
%     This function extracts all usefull time limits from whoAmI.m
%
%     INPUT
%        basedir     :  unit directory
%
%     OUTPUT   
%        TL          :  structure with all time limits
%
%
%   created by: 
%        Johannes Becherer
%        Tue Mar 27 12:40:09 PDT 2018
%


tl_default   = [datenum(2000, 1,1) datenum(2050, 1,1)];  % default tl


% gust of chipod
load([basedir '/calib/header.mat']);
is_chipod   =  isfield( head.coef, 'T1' );


%_____________________initialize______________________

TL.master   =  tl_default;
TL.P        =  tl_default;
TL.pitot    =  tl_default;
TL.acc      =  tl_default;
TL.cmp      =  tl_default;

if is_chipod
   TL.T1       =  tl_default;
   TL.Tp1      =  tl_default;
   TL.T2       =  tl_default;
   TL.Tp2      =  tl_default;
else
   TL.T        =  tl_default;
   TL.Tp       =  tl_default;
end


% check if whoAMIxist?
fname =  [basedir '/mfiles/whoAmI.m'];

%_____________________get from whoAmI______________________
if exist(fname)

 %---------------------read whoAmI----------------------
     startdir = pwd;

     % jump to whoAmI directory
     cd([basedir '/mfiles/']);
      
     DB = whoAmI;

     % jump back to inital directory
     cd(startdir);



 %---------------------convert whoAmI into time lims----------------------
  TL.master(1)       =  convert_timeStr2datenum( DB.instruments.start, tl_default(1));
  TL.master(2)       =  convert_timeStr2datenum( DB.instruments.stop, tl_default(2));

  tl_default = TL.master;


  %_____________________set all sub flags to master______________________
   TL.P        =  TL.master;
   TL.pitot    =  TL.master;
   TL.acc      =  TL.master;
   TL.cmp      =  TL.master;




  % chipod
   if is_chipod

      TL.T1       =  TL.master;
      TL.Tp1      =  TL.master;
      TL.T2       =  TL.master;
      TL.Tp2      =  TL.master;

      TL.T1(2)    =  convert_timeStr2datenum( DB.chipods.T1_stop, tl_default(2));
      TL.T2(2)    =  convert_timeStr2datenum( DB.chipods.T2_stop, tl_default(2));
      TL.Tp1(2)   =  convert_timeStr2datenum( DB.chipods.Tp1_stop, tl_default(2));
      TL.Tp2(2)   =  convert_timeStr2datenum( DB.chipods.Tp2_stop, tl_default(2));

      TL.P(2)     =  convert_timeStr2datenum( DB.chipods.P_stop, tl_default(2));
      TL.pitot(2) =  convert_timeStr2datenum( DB.chipods.pitot_stop, tl_default(2));
      TL.acc(2)   =  convert_timeStr2datenum( DB.chipods.acc_stop, tl_default(2));
      TL.cmp(2)   =  convert_timeStr2datenum( DB.chipods.cmp_stop, tl_default(2));

      % there's no way the differentiator returns good data if the T sensor is bad
      if TL.Tp1(2) > TL.T1(2)
          TL.Tp1(2) = TL.T1(2);
      end
      if TL.Tp2(2) > TL.T2(2)
          TL.Tp2(2) = TL.T2(2);
      end

   else  % gust

      TL.T       =  TL.master;
      TL.Tp      =  TL.master;

      TL.T(2)    =  convert_timeStr2datenum( DB.gust.T_stop, tl_default(2));
      TL.Tp(2)   =  convert_timeStr2datenum( DB.gust.Tp_stop, tl_default(2));

      TL.P(2)     =  convert_timeStr2datenum( DB.gust.P_stop, tl_default(2));
      TL.pitot(2) =  convert_timeStr2datenum( DB.gust.pitot_stop, tl_default(2));
      TL.acc(2)   =  convert_timeStr2datenum( DB.gust.acc_stop, tl_default(2));
      TL.cmp(2)   =  convert_timeStr2datenum( DB.gust.cmp_stop, tl_default(2));

      % there's no way the differentiator returns good data if the T sensor is bad
      if TL.Tp(2) > TL.T(2)
          TL.Tp(2) = TL.T(2);
      end

   end


end

end


function [t] = convert_timeStr2datenum( timeStr, defaultTime )

   if ~(isempty(timeStr))
         try
             try 
                t = datenum( timeStr);
             catch
                t = datenum(timeStr,'yyyy,mm,dd,HH,MM,SS');
             end
         catch
             t = defaultTime;
         end
         return;
   end

   t = defaultTime;
   

end

