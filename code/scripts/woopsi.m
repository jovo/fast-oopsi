% SI_ANALYSIS Test Josh Vogelstein's spike inference algorithms
% given a standard data structure. Plot electrophysiological data
% and calcium transients against with the inferred spike times from
% both the fast and slow spike inference techniques (oopsi).
%
% tamachado

clear; clear global; close all; clc

dataset = load('~/Research/oopsi/meta-oopsi/data/rafa/adam/2008/Imaging-SNR-Data.mat');
names = fieldnames(dataset);

% set simulation metadata
T               = 2930;     % # of time steps
V.dt            = 1/200;    % time step size
V.fast_iter_max = 1;        % # iterations of EM to estimate params
V.fast_plot     = 1;        % whether to plot with each iteration
V.save          = 1;        % whether to save results

% initialize params
P.a     = 1;                % scale
P.b     = 0;              % bias
tau     = 1;                % decay time constant for each cell
P.gam   = 1-V.dt/tau;       % set gam
P.lam   = 1;              % rate
P.sig   = 0.3;             % standard deviation

% simulate data
V.n = poissrnd(P.lam*V.dt*ones(T,1));
spt = find(V.n);
V.C = filter(1,[1 -P.gam],V.n);         % calcium concentration
V.F = P.a*(V.C+repmat(P.b,T,1))'+P.sig*randn(1,T);

%%
V.fast_iter_max = 10;
V.fast_plot     = 1;
V.est_sig       = 1;
V.name          = 'woopsi';                                 % filename for figure
V.name_dat      = ['../../data/' V.name];
V.name_fig      = ['../../figs/' V.name];

for i=8;
    cc      = dataset.(char(names(i)));
    F{i}    = cc.Fluorescence;
    V.dt    = median(diff(cc.FluorescenceTime));
    [n_f{i} P] = fast_oopsi(V.F,V); %  f.n;
    n_w{i}  = wiener_oopsi(V.F',V.dt,P); %  s.E.nbar;
    volt1{i} = cc.chanDev1_ai0_VoltageCh1;
    volt{i} = interp1(cc.time,volt1{i},cc.FluorescenceTime);
    temp    = GetSpikeTimes(volt1{i},0.7);
    n_t1{i}  = 0*volt1{i}; n_t1{i}(temp)=1;
    n_t{i}  = interp1(cc.time,n_t1{i},cc.FluorescenceTime);
    save(V.name_dat,'F','n_f','n_w','volt')
end

%%
j       = 8;
fig     = figure(2); clf,
V.T     = length(V.F);
nrows   = 4;
gray    = [.75 .75 .75];            % define gray color
inter   = 'tex';                    % interpreter for axis labels
xlims   = [100 930];                  % xmin and xmax for current plot
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
% plot(tvec_o,z1(F{j}(tvec_o)),'-k','LineWidth',lw);
plot(z1(V.F),'-k','LineWidth',lw);
ylab=ylabel([{'fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis('tight'); %[xlims 0 1.1])

% plot voltage data
i=i+1; h(i)=subplot(nrows,1,i); hold on
% plot(V.n,'-k','LineWidth',lw);
stem(spt,V.n(spt),'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'spike'}; {'train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([0 length(V.n) 0 1.1])

% plot fast-oopsi output
i=i+1; h(i)=subplot(nrows,1,i); hold on,
n_f{j}=n_f{j}/max(n_f{j});
spts=find(n_f{j}>1e-3);
stem(spts,n_f{j}(spts),'Marker','none','LineWidth',sw,'Color','k')
stem(spt,V.n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
axis([xlims 0 1.1])
hold off,
ylab=ylabel([{'fast'}; {'filter'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
box off

% plot wiener-oopsi output
i=i+1; h(i)=subplot(nrows,1,i); hold on,
wn=find(n_w{j}<0);
stem(wn,n_w{j}(wn)/max(n_w{j}),'Marker','none','LineWidth',sw,'Color',gray)

wp=find(n_w{j}>0);
stem(wp,n_w{j}(wp)/max(n_w{j}),'Marker','none','LineWidth',sw,'Color','k')

stem(spt,V.n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
axis([xlims min(n_w{j}/max(n_w{j})) 1.1])
hold off,
ylab=ylabel([{'wiener'}; {'filter'}],'Interpreter',inter,'FontSize',fs);
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