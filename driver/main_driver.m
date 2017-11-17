%     this is the driver for the main chi-processing
%
%    created:
%     Johannes Becherer
%     Tue Aug 23 12:26:41 PDT 2016

clear all;
close all;

%_____________________check processing flags or process?______________________
   do_dry_run = 0;

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%_____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);
   

%_____________________set processing flags______________________
   pflag = chi_processing_flags;     % get list of processing flags

   %---------------------gust or chipod----------------------
   if fids{1}(end-3)=='g' | fids{1}(end-3)=='G' % GusT
      pflag = pflag.c_gst(1);
   else                % chipod
      pflag = pflag.c_gst(0);
   end

   %---------------set processing flags automatically----------
   pflag = pflag.auto_set(basedir);

   %---------------------add manual flags----------------------
    %pflag = pflag.c_T1(0);       % switch off T1 if bad
    %pflag = pflag.c_T2(0);       % switch off T2 if bad

    %pflag = pflag.c_ic(1);       % switch on ic processing (default off)
    %pflag = pflag.c_vc(0);      % switch off viscous convective processing (default on)
    %pflag.master.epsp = 1;       % switch on eps calculation from pitot (default on)
  
    %pflag = pflag.c_vel_p(0);    % use pitot velocities 
    %pflag = pflag.c_vel_m(0);    % use mooring velocities 
    %pflag = pflag.c_Tzi(0);      % use local (interal) stratification 
    %pflag = pflag.c_Tzm(0);      % use mooring stratification 


    pflag = pflag.make_cons();   % make sub-flags consitent with master flags 

    %pflag.master.pumped = 1;    % do you use a surface pumped mooring (default = 0)

    pflag.master.use_compass = 1; % if 0, assume chipod vane moves it
                                  % into the flow perfectly

    pflag.master.use_pres = 0; % set to 0, accelerometers are
                               % working, else differentiate
                               % pressure to get speed past sensor
    pflag.master.parallel = 1;

   %---------------------get flag status----------------------
   pflag.status();



%_____________processing loop through all raw files__________________
if ~do_dry_run   
   % init parallel pool
   if(pflag.master.parallel)
      parpool;
      % parallel for-loop
      parfor f=1:length(fids)
            disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
         try % take care if script crashes that the parpoo is shut down
            chi_main_proc(basedir, fids{f}, pflag);
         catch
            disp(['!!!!!! ' fids{f} ' crashed while processing  !!!!!!' ]);
         end
      end
      % close parpool
      delete(gcp);
   else
      for f=1:length(fids)
         disp(['processing day ' num2str(f) ' of ' num2str(length(fids))]);
         chi_main_proc(basedir, fids{f}, pflag);
      end
   end


%_____________________merge all days______________________
disp('merge all days')
   %_loop through all processing flags for chi processing_
   for i = 1:length(pflag.id)
         [id, ~, ~, ~] = pflag.get_id(i);
         if pflag.proc.(id) % check if flag is active 
            disp([ id ' is being merged '  ]);
            ddir = ['chi' filesep 'chi_' id];
            % keep averaging window 0 here.
            % Only merge, average later in combine_turbulene.m
            chi_merge_and_avg(basedir, ddir, 0);
         end
   end

   % merge eps data
   if pflag.master.epsp
      disp('Pitot epsilon data are being merged')
      % keep averaging window 0 here.
      % Only merge, average later in combine_turbulene.m
      chi_merge_and_avg(basedir, 'eps', 0);
   end

   %_____________________combine all chi_data______________________
   % combine_turbulence;

end % ~dry run
