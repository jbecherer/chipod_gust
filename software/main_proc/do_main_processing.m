function [] = do_main_processing( basedir, pflag, time_lim )
%% [] = do_main_processing( basedir, [pflag], [time_lim] ) 
%  
%     This function drives the main processing for a given GusT or Chipod.
%
%     INPUT
%        basedir     :  base directory of the given instrument
%        pflag       :  processing flag object generated by chi_processing_flags.m
%                       if pflag is not given it will be generated automatically based on 
%                       inputs in the base directory
%        time_lim   :  time limits (notes this effects only the merged product not the processing time)
%
%   created by: 
%        Johannes Becherer
%        Thu Dec 28 12:33:57 PST 2017

ticstart = tic;
if nargin < 3 % if no time limits are set we use pratical no limits
   time_lim = [datenum(1900,1,1) datenum(2100,1,1)];
end

%____________________set automatic pflag______________________
   if nargin < 2

      pflag = chi_processing_flags;     % get list of processing flags

      %---------------set processing flags automatically----------
      pflag = pflag.auto_set(basedir);
      pflag.master.parallel = 0;
      %---------------------get flag status----------------------
      pflag.status();

   end

%_____________________do main processing______________________
   %_____________________get all raw files______________________

   [fids, fdate] = chi_find_rawfiles(basedir);


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
               ddir = ['chi' filesep 'chi_' id];
               % keep averaging window 0 here.
               % Only merge, average later in combine_turbulene.m
               chi_merge_and_avg(basedir, ddir, 0, time_lim );

               try
                   ddir = ['chi/chi_' id filesep 'stats'];
                   chi_merge_and_avg(basedir, ddir, 0);
               catch ME
                   disp(ME)
                   disp('Error! have the fitting stats been saved? Skipping...')
               end
            end
      end

   %_____________merge eps data______________________
   if pflag.master.epsp
      % keep averaging window 0 here.
      % Only merge, average later in combine_turbulene.m
      chi_merge_and_avg(basedir, 'eps', 0, time_lim);
   end

   disp('Finished running main processing.')
   toc(ticstart)
