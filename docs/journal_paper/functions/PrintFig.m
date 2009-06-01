function PrintFig(wh,DirName,FileName)

set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = [DirName FileName];
print('-depsc',FigName)
print('-dpdf',FigName)

end