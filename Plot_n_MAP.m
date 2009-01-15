function Plot_n_MAP(Pl,n)

% n=n(Pl.xlims(1):Pl.xlims(2));
% Pl.n=Pl.n(Pl.xlims(1):Pl.xlims(2));
% n=n/Pl.maxn;
% Pl.n=Pl.n/Pl.maxn;

cla, hold on
maxn(1) = max(n(Pl.xlims(1):Pl.xlims(2)));
maxn(2) = max(Pl.n(Pl.xlims(1):Pl.xlims(2)));

n=n/maxn(1);
Pl.n=Pl.n/maxn(2);

pos = find(n>0);                                % plot positive spikes in blue
stem(pos,n(pos),'Marker','none','LineWidth',Pl.sw,'Color',Pl.col(2,:))
neg = find(n<=0);                               % plot negative spikes in red
stem(neg,n(neg),'Marker','none','LineWidth',Pl.sw,'Color',Pl.col(1,:))

n(abs(n)<1e-3)=NaN;                             % set negligible values to NaN so they are not plotted
stem(Pl.n/maxn(2),'Marker','v','MarkerSize',Pl.vs,...   % plot real spike train
    'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');


axis([Pl.xlims min(min(n),0) 1])
ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:1,'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[])
set(gca,'XTickLabel',[])
box off
