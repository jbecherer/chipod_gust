function [fig] = plot_Turb(basedir, pflag, is_visible, runname)
% plot_Turb.m
%
% load Turb.m and plot all of the variables that have been processed

%%%%%%%%%%%%%%%%%%% DEFINE PARAMETERS %%%%%%%%%%%%%%%%%%%


% RUNNAME
% Name of folder in which Turb.mat is saved (folder within ../proc/)
% For example, if Turb.mat is in ../proc/Turb_min_dTdz_1e-3_60sec/Turb.mat,
% then runname = 'Turb_min_dTdz_1e-3_60sec'
% If Turb.mat is saved in ../proc/, leave runname = NaN;
if nargin < 4 
   runname = NaN;
end

if nargin < 3
    is_visible = 'on';
end

if nargin <2
% PROCESSING MODES TO PLOT
% There are up to 24 processing modes, which would you like to plot?
% Note: if any of these have not been processed, they will not plot
% automatically.
%_________ which estimates should I process?_______________________
   pflag = chi_processing_flags;

      %---------------------gust or chipod----------------------
      basedir = '../' 
      load([basedir '/calib/header.mat'])
      if isfield(head.coef, 'T') % GusT
         pflag = pflag.c_gst(1);
      else                % chipod
         pflag = pflag.c_gst(0);
      end

       pflag = pflag.auto_set(basedir);
      %---------------------add manual flags----------------------
       %pflag = pflag.c_T1(0);       % switch off T1 if bad
       %pflag = pflag.c_T2(0);       % switch off T2 if bad

       pflag = pflag.c_ic(1);       % switch on ic processing (default off)
       %pflag = pflag.c_vc(0);       % switch off viscous convective processing (default on)
       %pflag.master.epsp = 1;       % switch on eps calculation from pitot (default on)
     
       %pflag = pflag.c_vel_p(0);    % use pitot velocities 
       %pflag = pflag.c_vel_m(0);    % use mooring velocities 
       %pflag = pflag.c_Tzi(0);      % use local (interal) stratification 
       %pflag = pflag.c_Tzm(0);      % use mooring stratification 
      pflag = pflag.make_cons();     % make sub-flags consitent with master flags 
end

% set basedir 
if nargin <1
    here    =   pwd;                % mfiles folder
    basedir =   here(1:(end-6));    % substract the mfile folder
end

% load in Turb.mat
tic
if ~isnan(runname)
    disp(['loading ../proc/' runname '/Turb.mat'])
	load([basedir '/proc/' runname '/Turb.mat']);
else
    disp('loading ../proc/Turb.mat')
	load([basedir '/proc/Turb.mat']);
end
toc

ff = fields(Turb);
ff = {ff{1:end-1}}'; % remove readme structure

fig = CreateFigure(is_visible);
fig.Name = 'plot_Turb: compare all estimates';
[ax, ~] = create_axes(fig, 5, 1, 0);
      
 
%%%%%%%%%%%%%%%%% plot time series %%%%%%%%%%%%%%%%%
for a = 1:5
    % define parameters for each subplot
    if a == 1
        var = 'chi';
        labelstr = '\chi [K^2/s]';
        yl = [];
        yscale = 'log';
    elseif a == 2
        var = 'eps';
        labelstr = '\epsilon [m^2/s^3]';
        yl = [];
        yscale = 'log';
    elseif a == 3
        var = 'Kt';
        labelstr = 'K_t [m^2/s]';
        yl = [1e-7 1e0];
        yscale = 'log';
    elseif a == 4
        var = 'dTdz';
        labelstr = 'dTdz [C/m]';
        yl = [-1e-1, 0.3]; 
        yscale = 'linear';
    elseif a == 5
        var = 'spd';
        labelstr = '|u| [m/s]';
        yl = [0 2];
        yscale = 'linear';
    end
    
    % plot
    for f = 1:length(ff)
        if ~isstruct(Turb.(ff{f})) | strcmp( ff{f}, 'parameters' ) , continue; end
          if pflag.proc.(ff{f}) == 1
              if strcmpi(var, 'chi') | strcmpi(var, 'eps')
                  values = Turb.(ff{f}).(var);
                  Turb.(ff{f}).(var)(values == 0) = nan;
              end
            pj = f; p(pj) = plot(ax(a), Turb.(ff{f}).time, Turb.(ff{f}).(var),...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'));

            if isempty(yl)
                yl = [prctile(Turb.(ff{f}).(var), 1), nanmax(Turb.(ff{f}).(var))];
            end
        end
    end
    t = text_corner(ax(a), labelstr, 1);
    set(ax(a), 'Yscale', yscale);
    ylim(ax(a), sort(yl));
    set(ax(a),'tickdir','out')
    datetick(ax(a), 'keeplimits');
    
    % add some extras
    if a <4
        set(ax(a),'YTick',10.^(-10:1:0))
    end
    if a == 4
        set(ax(a), 'XTickLabels', [])
        plot(ax(a), xlim ,[0 0], 'k-')
        set(gcf, 'currentaxes', ax(a))
    end
    if a ~= 5
        set(ax(a),'XTickLabel',{' '})
    end

end

maybe_symlog(ax(4), 'y', -3);

linkaxes(ax, 'x');

%%%%%%%%%%%%%%%%% histograms %%%%%%%%%%%%%%%%%

ms = 8; %marker size of means and medians

squeeze_axes(ax, .8, 1)

pos = get(ax(1), 'Position');
axh(1) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.04 pos(4)] );
      
pos = get(ax(2), 'Position');
axh(2) = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.04 pos(4)] );

a=1;
hold(axh(a),'all');
yl = log10(get(ax(a), 'Ylim'));
bins = yl(1):diff(yl)/100:yl(2);
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1 | strcmp( ff{f}, 'parameters' ) , continue; end
      if pflag.proc.(ff{f}) == 1
        [Nchi,~] = histcounts( log10(Turb.(ff{f}).chi) , bins);
        pj = f; p(pj) = plot(axh(a), Nchi/nansum(Nchi) , bins(1:end-1)+diff(bins(1:2)*.5),...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'));
    end
end
xlims = get(axh(a),'XLim');
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1 | strcmp( ff{f}, 'parameters' ) , continue; end
       if pflag.proc.(ff{f}) == 1 & isfield(Turb.(ff{f}), 'stats')
        plot(axh(a), xlims(2)*0.9 ,log10( Turb.(ff{f}).stats.chimean),'<',...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'),'markersize',ms);
        plot(axh(a), xlims(2)*0.95 ,log10( Turb.(ff{f}).stats.chimedian),'+',...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'),'markersize',ms);
    end
end
ylim(axh(a), sort(yl));
set(axh(a), 'Xticklabel', {},'tickdir','out','YAxisLocation','right',...
    'YTick',-10:1:0)
xl = get(axh(a),'XLim');
text(axh(a),xl(2)*0.95,yl(2),' median -->','rotation',-90,'fontsize',10,...
    'horizontalalignment','left','verticalalignment','middle');
text(axh(a),xl(2)*0.9,yl(2),' mean -->','rotation',-90,'fontsize',10,...
    'horizontalalignment','left','verticalalignment','middle');
      
a=2;
hold(axh(a),'all');
yl = log10(get(ax(a), 'Ylim'));
bins = yl(1):diff(yl)/100:yl(2);
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1 | strcmp( ff{f}, 'parameters' ) , continue; end
       if pflag.proc.(ff{f}) == 1
        [Nchi,~] = histcounts( log10(Turb.(ff{f}).eps) , bins);
        pj = f; p(pj) = plot(axh(a), Nchi/nansum(Nchi) , bins(1:end-1)+diff(bins(1:2)*.5),...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'));
    end
end
xlims = get(axh(a),'XLim');
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1 | strcmp( ff{f}, 'parameters' ) , continue; end
       if pflag.proc.(ff{f}) == 1  & isfield(Turb.(ff{f}), 'stats')
        plot(axh(a), xlims(2)*0.9 , log10(Turb.(ff{f}).stats.epsmean),'<',...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'),'markersize',ms);
        plot(axh(a), xlims(2)*0.95 , log10(Turb.(ff{f}).stats.epsmedian),'+',...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'),'markersize',ms);
    end
end
ylim(axh(a), sort(yl));
% set(axh, 'Yticklabel', {}, 'Xticklabel', {},'tickdir','out')
set(axh(a), 'Xticklabel', {},'tickdir','out','YAxisLocation','right',...
    'YTick',-10:1:0)



%%%%%%%%%%%%%%%%% legend %%%%%%%%%%%%%%%%%
pos = get(ax(3), 'Position');
axl = axes('Position', [ pos(1)+pos(3) pos(2) 1-pos(3)-pos(1)-.04 pos(4)] );
hold(axl,'on');

clear nleg p ffstr
nleg = 1;
for f = 1:length(ff)
    if ~isstruct(Turb.(ff{f})) == 1 | strcmp( ff{f}, 'parameters' ) , continue; end
      if pflag.proc.(ff{f}) == 1
        pj = nleg; p(pj) = plot(axl, [0 1] ,[0 1],...
                'color', choose_color(ff{f},'color'), ...
                'LineWidth', choose_color(ff{f},'width'), ...
                'LineStyle', choose_color(ff{f},'style'));
        ffstr{nleg} = fix_underscore(ff{f});
        nleg = nleg+1;
    end
end
legend(p, ffstr);
set(axl, 'visible', 'off')
ylim(axl,[-1 -.5])





   
