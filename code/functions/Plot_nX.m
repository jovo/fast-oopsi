function Plot_nX(Pl,X)
hold on
for i=1:Pl.Nc
    plot(X(Pl.xlims(1):Pl.xlims(2),i),'Color',Pl.colors(i,:),'LineWidth',Pl.lw);
end
ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','Interpreter',Pl.interp,'FontName',Pl.font)
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims-Pl.xlims(1) min(X(:)) 1.1*max(X(:))])
box off
