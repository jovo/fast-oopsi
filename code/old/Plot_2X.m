function Plot_2X(Pl,X)
X1=X(Pl.xlims(1):Pl.xlims(2),1);
X2=X(Pl.xlims(1):Pl.xlims(2),2);
hold on
plot(X1,'Color',Pl.color1,'LineWidth',Pl.lw);
plot(X2,'Color',Pl.color2,'LineWidth',Pl.lw);
ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims-Pl.xlims(1) min(X(:)) max(X(:))])
box off
