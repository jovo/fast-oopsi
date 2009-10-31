% SI_ANALYSIS Test Josh Vogelstein's spike inference algorithms
% given a standard data structure. Plot electrophysiological data
% and calcium transients against with the inferred spike times from
% both the fast and slow spike inference techniques (oopsi).
%
% tamachado

clear; clear global; close all; clc

dataset = load('~/Research/oopsi/meta-oopsi/data/rafa/adam/2008/Imaging-SNR-Data.mat');
load '~/Research/oopsi/meta-oopsi/docs/posters/SFN09/oopsi.mat';
names = fieldnames(dataset);
%%
V.fast_do       = 1; 
V.smc_do        = 1; 
V.save          = 0; 
V.plot          = 1;
V.fast_iter_max = 10;
V.smc_iter_max  = 10;
V.N             = 99;
P.gam           = 1-(1/15)/1;
P.lam           = 10;

for i=8;
    cc      = dataset.(char(names(i)));
    F{i}    = cc.Fluorescence;
    V.dt    = median(diff(cc.FluorescenceTime));
%     [f s]   = run_oopsi(F{i},V,P);
    n_f{i}  = s.fastInferredSpikes{i}; %  f.n;
    n_s{i}  = s.smcInferredSpikes{i}; %  s.E.nbar;
    volt1{i} = cc.chanDev1_ai0_VoltageCh1;
    volt{i} = interp1(cc.time,volt1{i},cc.FluorescenceTime);
    temp    = GetSpikeTimes(volt1{i},0.7);
    n_t1{i}  = 0*volt1{i}; n_t1{i}(temp)=1;
    n_t{i}  = interp1(cc.time,n_t1{i},cc.FluorescenceTime);
    save('total','F','n_f','n_s','volt')
end


j       = i;
V.name_fig = '../../figs/smc_init';                                 % filename for figure
fig     = figure(1); clf,
V.T     = length(F{j});
nrows   = 4;
gray    = [.75 .75 .75];            % define gray color
inter   = 'tex';                    % interpreter for axis labels
xlims   = [70 V.T-15];              % xmin and xmax for current plot
fs      = 14;                       % font size
ms      = 5;                        % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
xticks  = 0:1/V.dt:V.T;             % XTick positions
skip    = round(length(xticks)/5);
xticks  = xticks(1:skip:end);
tvec_o  = xlims(1):xlims(2);        % only plot stuff within xlims

% plot fluorescence data
i=1; h(i)=subplot(nrows,1,i); hold on
plot(tvec_o,z1(F{j}(tvec_o)),'-k','LineWidth',lw);
ylab=ylabel([{'fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot voltage data
i=i+1; h(i)=subplot(nrows,1,i); hold on
plot(z1(volt{j}),'-k','LineWidth',lw);
ylab=ylabel([{'voltage'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([0 length(volt) 0 1.1])

% plot fast-oopsi output
i=i+1; h(i)=subplot(nrows,1,i); hold on,
n_f{j}=n_f{j}/max(n_f{j});
spts=find(n_f{j}>1e-3);
stem(spts,n_f{j}(spts),'Marker','none','LineWidth',sw,'Color','k')
stem(spts,n_t{j}(spts)+1,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
axis([xlims 0 1.1])
hold off,
ylab=ylabel([{'fast'}; {'filter'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
box off

% plot smc-oopsi output
i=i+1; h(i)=subplot(nrows,1,i); hold on,
spt=find(n_s{j}>1e-3);
stem(spt,n_s{j}(spt),'Marker','none','LineWidth',sw,'Color','k')
stem(spts,n_t{j}(spts)+1,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
axis([xlims 0 1.1])
hold off,
ylab=ylabel([{'smc'}; {'filter'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
box off

% label last subplot
set(gca,'XTick',xticks,'XTickLabel',round(xticks*V.dt*100)/100)
xlabel('Time (sec)','FontSize',fs)
linkaxes(h,'x')

% print fig
wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-depsc',V.name_fig)
print('-dpdf',V.name_fig)
saveas(fig,V.name_fig)