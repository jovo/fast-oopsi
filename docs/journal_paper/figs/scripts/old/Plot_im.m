function Plot_im(Pl,I)

imagesc(reshape(I,Pl.w,Pl.h)),
ylab=ylabel(Pl.ylab,'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'xTick',[],'YTick',[])
tit=title(Pl.tit);
set(tit,'FontSize',Pl.fs)