function Plot_X(Pl,X)
X=X(Pl.xlims(1):Pl.xlims(2));
plot(X,'Color',Pl.color,'LineWidth',Pl.lw);
ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims-Pl.xlims(1) min(X) max(X)])
box off
