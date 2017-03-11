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


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

%_____________________for automated PMEL mooring processing____________
    if use_pmel
        pmeldir = '~/TaoTritonPirataRama/'; % directory with pmel mooring files
                                            % (can obtain an updated copy from ganges)
        % which high-freq data file should I use?
        % 2m/10m/30m/hr
        velfreq = '30m';
        Tfreq = '10m';
        Sfreq = 'dy';

        % find start and end of depoyment from raw files
        rawdir = [basedir filesep 'raw' filesep];
        data = raw_load_chipod([rawdir fids{1}]);
        deployStart = data.datenum(1);
        data = raw_load_chipod([rawdir fids{end}]);
        deployEnd = data.datenum(end);

        % chipod location (positive North, East & Down)
        ChipodLon = 90; ChipodLat = 12; ChipodDepth = 15;
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
                          T2.time, T2.z, T2.T, T2.S, sdir);

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
