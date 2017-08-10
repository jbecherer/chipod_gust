do_parallel = 1;

%_____________________include path of processing files______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


% determine if chipod or gusT
load ../calib/header.mat
if isfield(head.coef, 'T1')
    isChipod = 1;
else
    isChipod = 0;
end


%____________________set directories______________________
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name

%_____________________get list of all raw data______________________
   [fids, fdate] = chi_find_rawfiles(basedir);

%_____________________average all Pitot raw data to 10 min bins______________________
disp('averaging Pitot raw  data in ./proc/Praw.mat ')
% init parallel pool
if(do_parallel)
    parpool;
    % parallel for-loop
    parfor f=1:length(fids)
        try % take care if script crashes that the parpoo is shut down
            disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);

            [~] = pitot_avg_raw_data(basedir, fids{f}, isChipod);
        catch ME
            disp(['!!!!!! ' fids{f} ' crashed while processing vel_p structure !!!!!!' ]);
            disp(['f = ' num2str(f)])
            disp(ME)
        end
    end
    % close parpool
    delete(gcp);
else
    for f=1:length(fids)
        disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
        [Praw] = pitot_avg_raw_data(basedir, fids{f}, isChipod);
    end
end

%_____________________merge individual files______________________
chi_merge_and_avg(basedir, 'Praw', 0);
