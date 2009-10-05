function Plot_2n(Pl,n)

cla, hold on
stem(n(:,1),'Marker','none','LineWidth',Pl.sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k');
if Pl.Nc>1
    stem(-n(:,2),'Marker','none','LineWidth',Pl.sw,'Color',Pl.gray,'MarkerFaceColor',Pl.gray,'MarkerEdgeColor',Pl.gray);
end
ylab=ylabel([{'Spike'}; {'Train'}],'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims min(-n(Pl.xlims(1):Pl.xlims(2),end)) max(n(Pl.xlims(1):Pl.xlims(2),1))])
% axis([Pl.xlims 0 max(n(Pl.xlims(1):Pl.xlims(2)))])
box off