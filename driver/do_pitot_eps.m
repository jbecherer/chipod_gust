%% this script does all the pitot epsilon processing
%     in ./proc/pitot_eps.mat
%  
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:01:26 PDT 2017


do_parallel = 0;     % use paralelle computing 
do_raw_proc = 0;     % do the raw data processing 
   save_spec = 0;    % shall the spectrum be saved
do_plot     = 1;     % generate a over view plot 
do_merge    = 0;     % merge provcessed files 


%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location

% get time limits from whoAmI;
   [TL] =   whoAmI_timeLimits(basedir);
   time_range      = TL.pitot;
   %time_range      = [datenum(2017,09,20) datenum(2017,09,21)];


if do_raw_proc
    spec_length   = 5; % in seconds 
    accfilter     = 1; % shall the spectrum be filtered Correlation with Acc
                       % as a result the spectrum will be shorter (needs longer spec_length)
                       % and processing time longer

    frange        =  [1/spec_length 2];


    generate_pitot_eps(basedir, do_parallel,  time_range, spec_length, frange,  save_spec, accfilter);
end

if do_merge
   disp('---- merging processed days ----')
   chi_merge_and_avg(basedir, 'pitot_eps', 0, time_range);
end


if do_plot
      
      if save_spec
         peps_spec_plot;
      end

      if isfield(TL, 'T') % gust
         whichTurb = 'pmg';
      else  % chipod
         if TL.Tp1(2) < TL.Tp2(2)
            whichTurb = 'pm2';
         else
            whichTurb = 'pm1';
         end
      end
         
      if exist([basedir '/proc/Turb.mat'])
         [fig, fig1] = plot_compare_PepsTurb( basedir, whichTurb)
            print(fig,['../pics/Pitot_eps_Turb.png'],'-dpng','-r200','-painters')
            print(fig1,['../pics/Pitot_eps_Turb_hist2d.png'],'-dpng','-r200','-painters')
      end

      [fig] = plot_pitot_eps( basedir)
      print(fig,['../pics/Pitot_eps_bbl.png'],'-dpng','-r200','-painters')


            
end
