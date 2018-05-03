function [TL] =   db_timeLimits(db_index, p2g)
%%    [TL] =   db_timeLimits(db_index, p2g)
%
%     This function extracts all usefull time limits from omg.sqlite
%
%     INPUT
%        db_index     : index of instrument in the instrument table
%        p2g          : path 2 ganges
%
%     OUTPUT   
%        TL          :  structure with all time limits
%
%
%   created by: 
%        Johannes Becherer
%        Tue Mar 27 12:40:09 PDT 2018
%

if nargin < 2
   p2g = get_ganges_path();
end
path2database = [p2g '/work/database/omg.sqlite'];


tl_default   = [datenum(2000, 1,1) datenum(2050, 1,1)];  % default tl

%_____________________initialize______________________

TL.master   =  tl_default;
TL.P        =  tl_default;
TL.pitot    =  tl_default;
TL.acc      =  tl_default;
TL.cmp      =  tl_default;

% chipods
TL.T1       =  tl_default;
TL.Tp1      =  tl_default;
TL.T2       =  tl_default;
TL.Tp2      =  tl_default;
% gusts
TL.T        =  tl_default;
TL.Tp       =  tl_default;


%_____________________get info from instrument table____________________

     % open data base
     omg_db = sqlite(path2database);
     sqlstr =  ['select name, type, start, stop from instruments where id = ' num2str(db_index)];
     try
         inst_data = fetch(omg_db, sqlstr );
     end
     % close database
     close(omg_db);

     DB.instruments.name    =  char(inst_data{1,1});
     DB.instruments.type    =  char(inst_data{1,2});
     DB.instruments.start   =  char(inst_data{1,3});
     DB.instruments.stop    =  char(inst_data{1,4});


% gust of chipod
if contains(DB.instruments.type, 'chipod' ,'IgnoreCase',true) 
   is_chipod   =  1;
elseif contains(DB.instruments.type, 'gust' ,'IgnoreCase',true) 
   is_chipod   =  0;
else
   warning(['The db_index corresponds to ' DB.instruments.name ', which is not a chipod or gust.'])
   disp('The time limits wil be set to default!')
   return
end

%_____________________initialize______________________

if is_chipod
   TL	=	rmfield( TL, 'T');
   TL	=	rmfield( TL, 'Tp');
else
   TL	=	rmfield( TL, 'T1');
   TL	=	rmfield( TL, 'Tp1');
   TL	=	rmfield( TL, 'T2');
   TL	=	rmfield( TL, 'Tp2');
end




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

		% read chipod table		
		% open data base
		omg_db = sqlite(path2database);
		sqlstr =  ['select T1_stop, T2_stop, Tp1_stop, Tp2_stop, P_stop, ' ...
		 			 'pitot_stop, acc_stop, cmp_stop from chipods where inst_id = ' num2str(db_index)];
		try
		 	chipod_data = fetch(omg_db, sqlstr );
		end
		% close database
		close(omg_db);

		DB.chipods.T1_stop		=  char(chipod_data{1,1});
		DB.chipods.T2_stop      =  char(chipod_data{1,2});
		DB.chipods.Tp1_stop     =  char(chipod_data{1,3});
		DB.chipods.Tp2_stop     =  char(chipod_data{1,4});
		DB.chipods.P_stop			=  char(chipod_data{1,5});
		DB.chipods.pitot_stop   =  char(chipod_data{1,6});
		DB.chipods.acc_stop     =  char(chipod_data{1,7});
		DB.chipods.cmp_stop     =  char(chipod_data{1,8});


      TL.T1(2)    =  convert_timeStr2datenum( DB.chipods.T1_stop, tl_default(2));
      TL.T2(2)    =  convert_timeStr2datenum( DB.chipods.T2_stop, tl_default(2));
      TL.Tp1(2)   =  convert_timeStr2datenum( DB.chipods.Tp1_stop, tl_default(2));
      TL.Tp2(2)   =  convert_timeStr2datenum( DB.chipods.Tp2_stop, tl_default(2));

      TL.P(2)     =  convert_timeStr2datenum( DB.chipods.P_stop, tl_default(2));
      TL.pitot(2) =  convert_timeStr2datenum( DB.chipods.pitot_stop, tl_default(2));
      TL.acc(2)   =  convert_timeStr2datenum( DB.chipods.acc_stop, tl_default(2));
      TL.cmp(2)   =  convert_timeStr2datenum( DB.chipods.cmp_stop, tl_default(2));

      % there's no way the differentiator gives good data once the T sensor dies
      if TL.Tp1(2) > TL.T1(2)
          TL.Tp1(2) = TL.T1(2);
      end
      if TL.Tp2(2) > TL.T2(2)
          TL.Tp2(2) = TL.T2(2);
      end

   else  % gust

      TL.T       =  TL.master;
      TL.Tp      =  TL.master;

		% read gust table		
		% open data base
		omg_db = sqlite(path2database);
		sqlstr =  ['select T_stop, Tp_stop, P_stop, ' ...
		 			 'pitot_stop, acc_stop, cmp_stop from gust where inst_id = ' num2str(db_index)];
		try
		 	gust_data = fetch(omg_db, sqlstr );
		end
		% close database
		close(omg_db);

		DB.gust.T_stop			=  char(gust_data{1,1});
		DB.gust.Tp_stop     	=  char(gust_data{1,2});
		DB.gust.P_stop			=  char(gust_data{1,3});
		DB.gust.pitot_stop   =  char(gust_data{1,4});
		DB.gust.acc_stop     =  char(gust_data{1,5});
		DB.gust.cmp_stop     =  char(gust_data{1,6});

      TL.T(2)    =  convert_timeStr2datenum( DB.gust.T_stop, tl_default(2));
      TL.Tp(2)   =  convert_timeStr2datenum( DB.gust.Tp_stop, tl_default(2));

      TL.P(2)     =  convert_timeStr2datenum( DB.gust.P_stop, tl_default(2));
      TL.pitot(2) =  convert_timeStr2datenum( DB.gust.pitot_stop, tl_default(2));
      TL.acc(2)   =  convert_timeStr2datenum( DB.gust.acc_stop, tl_default(2));
      TL.cmp(2)   =  convert_timeStr2datenum( DB.gust.cmp_stop, tl_default(2));

      if TL.Tp(2) > TL.T(2)
          TL.Tp(2) = TL.T(2);
      end

   end




end


function [t] = convert_timeStr2datenum( timeStr, defaultTime )

   if ~(isempty(timeStr))
         try 
            t = datenum( timeStr);
         catch
            t = defaultTime;
         end
         return;
   end

   t = defaultTime;
   

end

