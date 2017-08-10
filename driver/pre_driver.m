%     this is the driver for the prep-processing for  chi-processing
%
%   created by: 
%        Johannes Becherer
%        Tue Sep 20 10:51:19 PDT 2016

clear all;
close all;


%_____________________set processing flags______________________
   do_parallel = 0;     % use paralelle computing 
   do_temp     = 0;     % generate temp.mat 
   do_vel_p    = 0;     % generate vel_p.mat 
   do_vel_m    = 0;     % generate vel_m.mat
   do_dTdz_m   = 0;     % generate dTdz_m.mat
   do_dTdz_i   = 0;     % generate dTdz_i.mat 
   use_pmel    = 0;     % use TAO/TRITON/PIRATA/RAMA mooring data?
   use_mooring_sal = 0; % use mooring salinity along with dTdz_i
                        % to estimate N^2 in dTdz_i.
                        % otherwise code assumes fixed salinity=35.
   use_TS_relation = 0; % fit TS relation to estimate N2 from
                        % mooring data? Use (with caution) when you
                        % have only 1 salinity sensor
   modify_header = 0;   % if 1, specify header corrections below

   % declination - get values from https://www.ngdc.noaa.gov/geomag-web/#declination
   DeployDecl = 0; % at deployment location
   CorvallisDecl = 15+44/60; % at corvallis

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

%_____________________for automated PMEL mooring processing____________
    if use_pmel
        pmeldir = '~/ganges/data/TaoTritonPirataRama/'; % directory with pmel mooring files
                                            % (can obtain an updated copy from ganges)
        % which high-freq data file should I use?
        % 2m/10m/30m/hr
        velfreq   = '30m';
        Tfreq     = '10m';
        Sfreq     = 'dy';

        % find start and end of depoyment from raw files
        data         = raw_load_chipod([rawdir fids{1}]);
        deployStart  = data.datenum(1);
        data         = raw_load_chipod([rawdir fids{end}]);
        deployEnd    = data.datenum(end);

        % chipod location (positive North, East & Down)
        ChipodLon = 90; ChipodLat = 12; ChipodDepth = 15;
    end

%_______________create header files if necessary______________________

   % will create header.mat file if necessary
   % if header.mat exists, it will read it
   head = chi_get_calibration_coefs(basedir);

   % will create header.mat file if necessary
   % if header_p.mat exists, it will read it
   W  = chi_get_calibration_coefs_pitot(basedir);

%_____________________make header corrections if necessary______________

    if modify_header

        % _____usual calibrations______
        % head.coef.CMP(1) = 0;
        % save header in proper destination
        % fid = [basedir filesep 'calib' filesep 'header.mat'] ;
        % save(fid, 'head');

        % _____pitot calibrations______
        % W.T  = [0 -0.003154669 0 0 0];
        % W.Ps = [0 0 0 0 0]; % pressure sensor is bad; accounted for in offset
        % W.Tilt = [0 0.000088684 0 0 0];
        % W.Pd = [0 0.0003995 0 0 0]; %if slope>1 else W.Pd = [0 slope 0 0 0];
        % assert(W.Pd(2) < 1, 'WPd(2) > 1 !');
        % % offsets
        % W.V0 = 0;
        % W.P0 = 0;
        % W.T0 = 0;
        % save header in proper destination
        % fid = [basedir filesep 'calib' filesep 'header_p.mat'] ;
        % save(fid, 'W');

    end

%_____________________account for declination__________________________

    [~, headtest] = raw_load_chipod([basedir filesep 'raw' filesep fids{1}]);

    % if these are equal, declination correction hasn't been made yet
    if isequal(head.coef.CMP(1), headtest.coef.CMP(1))
        disp('accounting for declination ...');
        % chi_calibrate_chipod adds head.coef.CMP(1) to raw_data.CMP/10
        % hence, we need to change sign here.
        head.coef.CMP(1) = -head.coef.CMP(1) - CorvallisDecl + DeployDecl;
        fid = [basedir filesep 'calib' filesep 'header.mat'] ;
        save(fid, 'head');
    else
        if isequal(head.coef.CMP(1), -headtest.coef.CMP(1) - CorvallisDecl + DeployDecl)
            disp('declination has already been accounted for.')
        end
    end

%%%%%%%%%%%%%%%%%%% temp processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_temp
   %_____________processing loop through all raw files__________________

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_T_proc(basedir, fids{f});
            catch
               disp(['!!!!!! ' fids{f} ' crashed while processing T structure !!!!!!' ]);
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            chi_T_proc(basedir, fids{f});
         end
      end

   %____________________merge individual files______________________
      
      % average 20 sec
      chi_merge_and_avg(basedir, 'temp', 20);
end


%%%%%%%%%%%%%%%%%%% generating Pitot velocity input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_vel_p

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
      disp([fid ' does not exist. Run calibrate_pitot first !']);
   end

   vel_p.time  = P.time;
   vel_p.spd   = P.spd;
   vel_p.U     = P.U;
   vel_p.u     = real(P.U);
   vel_p.v     = imag(P.U);

   save('../input/vel_p.mat', 'vel_p');
   
end


%%%%%%%%%%%%%%%%%%% mooring velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_vel_m
    sdir  = [basedir filesep 'input' filesep];

    if use_pmel
        moor = ExtractUVFromTaoTritonPirataRama(ChipodLon, ChipodLat, ...
                                                ChipodDepth, deployStart, ...
                                                deployEnd, pmeldir, ...
                                                'RAMA', velfreq);
    end

    %_______ EXAMPLE________________
    % load('../../../mooring_data/mooring_Pirata14_524.mat') ;

    chi_generate_vel_adcp(moor.time, moor.depth, moor.u, moor.v, moor.depth, sdir);
end


%%%%%%%%%%%%%%%%%%% mooring dTdz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_dTdz_m
      sdir  = [basedir filesep 'input' filesep];

      if use_pmel
          [T1, T2] = ExtractTSFromTaoTritonPirataRama(ChipodLon, ChipodLat, ...
                                                      ChipodDepth, deployStart, ...
                                                      deployEnd, pmeldir, 'RAMA', ...
                                                      Tfreq, Sfreq);
      end

 %_______ EXAMPLE________________
      %  load('../../G002/proc/temp.mat') ; % surounding instruments
      %     T1.time = T.time; 
      %     T1.z    = nanmedian(T.depth); 
      %     T1.T    = T.T; 
      %     T1.S    = ones(size((T.T)))*35; 
      %  load('../../G011/proc/temp.mat') ; % surounding instruments
      %     T2.time = T.time; 
      %     T2.z    = nanmedian(T.depth); 
      %     T2.T    = T.T; 
      %     T2.S    = ones(size((T.T)))*35; 

      chi_generate_dTdz_m(T1.time, T1.z, T1.T, T1.S, ...
                          T2.time, T2.z, T2.T, T2.S, sdir, use_TS_relation);

      %__________________recalculate N^2 using processed mooring salinity____________________

      if use_mooring_sal
          if ~exist('../input/dTdz_i.mat', 'file')
              error(['Create dTdz_i.mat first. Run pre_driver with ' ...
                     'do_dTdz_i=1']);
          end

          load ../input/dTdz_i.mat
          load ../input/dTdz_m.mat

          % interpolate to Tz_i.time
          dSdz = interp1(T1.time, (T1.S-T2.S)/abs(T1.z-T2.z), Tz_i.time);
          Smean = interp1(T1.time, (T1.S + T2.S)/2, Tz_i.time);

          Tnames = {'T1', 'T2', 'T12'};
          Tznames = {'Tz1', 'Tz2', 'Tz12'};
          Nnames = {'N2_1', 'N2_2', 'N2_12'};
          for ii=1:length(Tnames)
              Tii = Tz_i.(Tnames{ii});
              dTdz = Tz_i.(Tznames{ii});

              alpha = sw_alpha(Smean, Tii, ChipodDepth);
              beta = sw_beta(Smean, Tii, ChipodDepth);

              Tz_i.(Nnames{ii}) = -9.81 * (-alpha.*dTdz + beta.*dSdz);
          end

          Tz_i.S = Smean;
          Tz_i.Sz = dSdz;
          save('../input/dTdz_i.mat', 'Tz_i');
      end
end

%%%%%%%%%%%%%%%%%%% internal dTdz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
if do_dTdz_i
   %_____________processing loop through all raw files__________________

      dt   = 60; % sec bits of data for analysis
      do_P = 0; % use pressure instead of acceleration to get z 

      disp('calculating the intrenal dTdz');
      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_generate_dTdz_i(basedir, fids{f}, dt, do_P);
            catch
               disp(['!!!!!! ' fids{f} ' crashed while processing  internal dTdz structure !!!!!!' ]);
            end
         end
         % close parpool
         delete(gcp);
      else
         for f=1:length(fids)
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
            chi_generate_dTdz_i(basedir, fids{f}, dt, do_P);
         end
      end

   %____________________merge individual files______________________
      chi_merge_and_avg(basedir, 'dTdz', 600);

   %_____________________cp result to the input directory______________________
   ! cp ../proc/dTdz.mat ../input/dTdz_i.mat
end
