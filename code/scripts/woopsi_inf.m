clear; clear global; 

% set algorithm variables
V.T             = 700;     % # of time steps
V.dt            = 1/60;    % time step size
V.fast_iter_max = 1;        % # iterations of EM to estimate params
V.fast_plot     = 1;        % whether to plot with each iteration
V.save          = 1;        % whether to save results
V.fast_iter_max = 1;
V.fast_plot     = 0;
V.est_sig       = 1;
V.name          = 'woopsi_inf';                                 % filename for figure
V.name_dat      = ['../../data/' V.name];
V.name_fig      = ['../../figs/' V.name];

% initialize params
P.a     = 1;                % scale
P.b     = 0;              % bias
tau     = 1;                % decay time constant for each cell
P.gam   = 1-V.dt/tau;       % set gam
P.lam   = 1;              % rate
P.sig   = 0.3;             % standard deviation

% simulate data
V.n = poissrnd(P.lam*V.dt*ones(V.T,1));
spt = find(V.n);
V.C = filter(1,[1 -P.gam],V.n);         % calcium concentration
V.F = P.a*(V.C+repmat(P.b,V.T,1))'+P.sig*randn(1,V.T);

n_fast   = fast_oopsi(V.F,V,P);
n_wiener = wiener_oopsi(V.F,V.dt,P); %  s.E.nbar;

%%
fig     = figure(2); clf,
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
plot(z1(V.F),'-k','LineWidth',lw);
ylab=ylabel([{'fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis('tight'); %[xlims 0 1.1])

% plot voltage data
i=i+1; h(i)=subplot(nrows,1,i); hold on
stem(spt,V.n(spt),'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'spike'}; {'train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([0 length(V.n) 0 1.1])

% plot fast-oopsi output
i=i+1; h(i)=subplot(nrows,1,i); hold on,
n_fast=n_fast/max(n_fast);
spts=find(n_fast>1e-3);
stem(spts,n_fast(spts),'Marker','none','LineWidth',sw,'Color','k')
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
wn=find(n_wiener<0);
stem(wn,n_wiener(wn)/max(n_wiener),'Marker','none','LineWidth',sw,'Color',gray)

wp=find(n_wiener>0);
stem(wp,n_wiener(wp)/max(n_wiener),'Marker','none','LineWidth',sw,'Color','k')

stem(spt,V.n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
axis([xlims min(n_wiener/max(n_wiener)) 1.1])
hold off,
ylab=ylabel([{'Wiener'}; {'filter'}],'Interpreter',inter,'FontSize',fs);
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