%% this script does all the basic calibration excluding chi-processing and save 1s infomation
%     in ./proc/temp.mat
%  
%
%   created by: 
%        Johannes Becherer
%        Wed Aug 16 16:01:26 PDT 2017


do_parallel = 0;     % use paralelle computing 
do_raw_proc = 0;     % do the raw data processing 
do_plot     = 1;     % generate a over view plot 

time_range = [datenum(2000, 1, 1, 0, 0, 0) ...
              datenum(2060, 1, 1, 0, 0, 0)];

dtind = 600; % every 10 minutes, assuming 1 second estimates

%_____________________include path of processing flies______________________
addpath(genpath('./chipod_gust/software/'));% include  path to preocessing routines


%____________________set directories______________________    
   here    =   pwd;                % mfiles folder
   basedir =   here(1:(end-6));    % substract the mfile folder
   savedir =   [basedir 'proc/'];  % directory directory to save data
   unit    = chi_get_unit_name(basedir); % get unit name
   rawdir       = [basedir filesep 'raw' filesep]; % raw files location


if do_raw_proc
   %_____________________get list of all raw data______________________
      [fids, fdate] = chi_find_rawfiles(basedir);


   %_____________processing loop through all raw files__________________

   disp('calibrating all base quantities in ./proc/temp.mat ')

      % init parallel pool
      if(do_parallel)
         parpool;
         % parallel for-loop
         parfor f=1:length(fids)
            try % take care if script crashes that the parpoo is shut down
               disp(['calculating file ' num2str(f) ' of ' num2str(length(fids))]);
               chi_T_proc(basedir, fids{f});
            catch ME
               disp(['!!!!!! ' fids{f} ' crashed while processing T structure !!!!!!' ]);
               disp(ME)
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
      % average 1 sec
      chi_merge_and_avg(basedir, 'temp', 0);


   %____________________create motion.mat file______________________

   load([basedir 'proc' filesep 'temp.mat']);

   motion.chiangle = atan2(T.a_vel_y, T.a_vel_x) * 180/pi;
   motion.time = T.time;
   motion.cmp = T.cmp;
   motion.a_vel_x = T.a_vel_x;
   motion.a_vel_y = T.a_vel_y;
   motion.a_vel_z = T.a_vel_z;

   save([basedir 'proc' filesep 'motion.mat'], 'motion', '-v7.3');


end

%_____________________make multiple summary plots______________________
if do_plot

  if ~exist('T', 'var')
      load([basedir 'proc' filesep 'temp.mat']);
  end

   tind(1) = find_approx(T.time, time_range(1), 1);
   tind(2) = find_approx(T.time, time_range(2), 1);
   tl = T.time(tind);

   col = get(groot,'DefaultAxesColorOrder');
   

    fig = figure('Color',[1 1 1],'visible','on','Paperunits','centimeters',...
            'Papersize',[30 30],'PaperPosition',[0 -1 30 30]);


            [ax, ~] = create_axes(fig, 6, 1, 0);
            
            a = 1;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.depth(1:dtind:end), 'color', [0  0 0 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['P [dbar]'], 7);
               
            a = 2;
            if isfield(T, 'T1') 
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T1(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
               pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T2(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                  t1 = text_corner(ax(a), ['T1 [^\circ C]'], 1);
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['T2 [^\circ C]'], 3);
                  t2.Color = [col(2,:)];
            else
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T(1:dtind:end), 'color', [0  0 0 1], 'Linewidth', 1);
                  xlim(ax(a), tl);
                  t = text_corner(ax(a), ['T [^\circ C]'], 1);
            end
               xlim(ax(a), tl);
               
            
            a = 3;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.AX(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
            pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.AY(1:dtind:end) + 4, ...
                                 'color', [col(pj,:) 1], 'Linewidth', 1);
            pj = 3; p(pj) = plot(ax(a), T.time(1:dtind:end), T.AZ(1:dtind:end) + 9.81 - 4, ...
                                 'color', [col(pj,:) 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t1 = text_corner(ax(a), ['AX [m s^{-2}]'], 1);
               t1.Color = [col(1,:)];
               t2 = text_corner(ax(a), ['AY (offset) [m s^{-2}]'], 2);
               t2.Color = [col(2,:)];
               t3 = text_corner(ax(a), ['AZ (offset) [m s^{-2}]'], 3);
               t3.Color = [col(3,:)];
               ylim(ax(a), [-1, 1]*8);
            
            a = 4;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.cmp(1:dtind:end), '.', 'color', [0  0 0 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['compass [deg]'], 7);

            a = 5;
            pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.W(1:dtind:end), 'color', [0  0 0 1], 'Linewidth', 1);
               xlim(ax(a), tl);
               t = text_corner(ax(a), ['pitot [volt]'], 7);

            a = 6;
            if isfield(T, 'T1') 
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T1Pt(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
               pj = 2; p(pj) = plot(ax(a), T.time(1:dtind:end), T.T2Pt(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                  xlim(ax(a), tl);
                  t1 = text_corner(ax(a), ['T1P [volt]'], 1);  
                  t1.Color = [col(1,:)];
                  t2 = text_corner(ax(a), ['T2P [volt]'], 3);  
                  t2.Color = [col(2,:)];
            else
               pj = 1; p(pj) = plot(ax(a), T.time(1:dtind:end), T.TPt(1:dtind:end), 'color', [col(pj,:) 1], 'Linewidth', 1);
                  xlim(ax(a), tl);
                  t1 = text_corner(ax(a), ['TP [volt]'], 1);  
                  t1.Color = [col(1,:)];
            end
               
            linkaxes(ax, 'x');   

            abc='abcdefghijklmnopqrst';
            for a = 1:(size(ax,1)*size(ax,2))
               text_corner(ax(a), abc(a), 9);
               set(ax(a), 'Xtick', ceil(tl(1)):round(diff(tl)/5):floor(tl(2)));
            end
            


            datetick(ax(a), 'keepticks',  'keeplimits');
            
            unit    = chi_get_unit_name(basedir); % get unit name
            t = text_corner(ax(1), ['unit ' unit], -2);

	%_____________________save pic______________________
	print(gcf,[basedir 'pics' filesep 'temp.png' ],'-dpng','-r200','-painters')
%	savefig(gcf,[basedir 'pics' filesep 'temp.fig' ])

        %% displacement histogram

        CreateFigure; hold on;
        histogram(sqrt(2)*movstd(T.a_dis_x(tind(1):tind(2)), 60), ...
                  'displaystyle', 'stairs', 'linewidth', 1.5)
        histogram(sqrt(2)*movstd(T.a_dis_y(tind(1):tind(2)), 60), ...
                  'displaystyle', 'stairs', 'linewidth', 1.5)
        histogram(sqrt(2)*movstd(T.a_dis_z(tind(1):tind(2)), 60), ...
                  'displaystyle', 'stairs', 'linewidth', 1.5)
        legend('x', 'y', 'z');
        xlabel('Displacement = sqrt(2) x std(integrated accel_{x,y,z}) over 1 minute')
        print(gcf,[basedir 'pics' filesep 'disp.png' ],'-dpng','-r200','-painters')
end
