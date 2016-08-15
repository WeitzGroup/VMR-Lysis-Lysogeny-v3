% Part of the code used in:
% Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
% 
% From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v2
% MIT License

clf;
clear all
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figlv_ratio_label';

tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = tmpfilename;
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);
%set(gcf,'Position', [221 434 1050 372]);
set(gcf,'Position',[343 321 586 484]);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');

fp=fopen('lv_stats.mat');
if (fp<=0)
  % main data goes here
  % Assigns parameters & variables
  clear info
  clear stats
  info.r=1/24; %hr^-1
  info.d=1/48; %hr^-1
  info.K=10^7; %cells/ml
  info.phi=10^-8;   %ml/cells/hr
  info.m=1/6;  %hr^-1
  info.beta=20;% burst size
  info0=info;
  x0_array=cell2mat(struct2cell(info));
  xl_array=x0_array/10;
  xu_array=x0_array*10;
    
  % Sample a large number of times
  % Do this once and save
  
  %Number of Resamples when doing Latin hypercube
  nruns = 10;
  %Number of Samples in Latin Hypercube
  nS = 10^3;
  
  count=0;
  more off
  for j=1:nruns,
    j
    %Sample between upper and lower bounds (uniform prob. distribution)
    %sampling uniformly in log-space
    LHSample = LHSmid(nS,xl_array,xu_array);
    %Loop through sampled points
    for k=1:size(LHSample,1)
        info = array2vstruct(LHSample(k,:),fieldnames(info0));
        theory.N=info.m/(info.beta*info.phi);
        theory.V=(info.r*(1-theory.N/info.K)-info.d)/info.phi;
        if (theory.N>0 & theory.V>0)
          count=count+1;
   	  stats.N(count)=theory.N;
          stats.V(count)=theory.V;
        end
    end
  end
  stats.LH=LHSample;
else
  load('lv_stats');
end

set(gca,'fontsize',20);
x=log10(stats.N);
y=log10(stats.V./stats.N);
[p,s]=polyfit(x,y,1);
tmph=loglog(10.^x,10.^p(2)*(10.^x).^p(1),'r-');
set(tmph,'linewidth',3,'color',[0.75 0.75 0.75]);
hold on
loglog(stats.N,stats.V./stats.N,'k.');
hold on
set(gca,'xtick',10.^[-1:1:10]);
set(gca,'ytick',10.^[-5:1:5]);
set(gca,'fontsize',20);
xlabel('Microbial cell density, $N^{\ast}$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('VMR, $V^{\ast}/N^{\ast}$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
hold on
mYX = [y',x'];
nxbins=25;
nybins=25;
minx = max(min(x),-2);
maxx = max(x);
miny = max(min(y),-2);
maxy = max(y);
stepx = (maxx-minx)/(nxbins+1);
stepy = (maxy-miny)/(nybins+1);
xc = (linspace(minx,maxx,nxbins+1))';
yc = (linspace(miny,maxy,nybins+1))';
nXBins = length(xc);
nYBins = length(yc);
vXLabel = 0.5*(xc(1:(nXBins-1))+xc(2:nXBins));
vYLabel = 0.5*(yc(1:(nYBins-1))+yc(2:nYBins));
mHist2d = hist2d_mat(mYX,yc,xc);
contour(10.^vXLabel,10.^vYLabel,mHist2d/sum(sum(mHist2d)),20);
set(gca,'fontsize',20);
colorbar
hold on
% Redraw again to have it lay on top of points
[p,s]=polyfit(x,y,1);
tmph=loglog(10.^x,10.^p(2)*(10.^x).^p(1),'r-');
set(tmph,'linewidth',3,'color',[0.75 0.75 0.75]);
tmph=loglog(10.^[2.5 8.5],[1 1],'r--');
set(tmph,'linewidth',3);
tmph=loglog(10.^[2.5 8.5],[100 100],'r--');
set(tmph,'linewidth',3);
tmpt=text(10^6.5,10^2.8,'Lotka-Volterra');
set(tmpt,'fontsize',14,'interpreter','latex');
tmpt=text(10^2.7,2,'1:1');
set(tmpt,'fontsize',14,'interpreter','latex');
tmpt=text(10^2.6,200,'100:1');
set(tmpt,'fontsize',14,'interpreter','latex');
xlim([10^2.5 10^8.5]);
ylim([10^-2 10^3]);


% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
%set(gca,'ytick',10.^[-2:2:10]);

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
tmps=sprintf('$VMR\\propto N^{%4.2f}$',p(1));
tmplh = legend(tmps,'location','SouthWest');
set(tmplh,'interpreter','latex','fontsize',16);
% remove box
%set(tmplh,'visible','off')
%set(tmplh,'legend','boxoff');
%legend('boxoff');
set(tmplh,'box','off');

% title('','fontsize',24)
% 'horizontalalignment','left');

% for writing over the top
% coordinates are normalized again to (0,1.0)
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
% first two points are normalized x, y positions
% text(,,'','Fontsize',14);

% automatic creation of postscript
% without name/date
psprintc(tmpfilenoname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*
more on
