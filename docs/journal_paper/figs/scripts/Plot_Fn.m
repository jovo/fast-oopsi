function Plot_Fn(F,P,Pl,i)

F_proj=P.a(:,i)\(F-repmat(P.b',Pl.T,1))';
plot(z1(F_proj)+1,'Color',Pl.c(i,:),'LineWidth',Pl.lw)
bar(Pl.n(:,i),'FaceColor',Pl.c(i,:),'EdgeColor',Pl.c(i,:),'LineWidth',2)
% stem(Pl.n(:,i),'Marker','v','LineStyle','none','MarkerFaceColor',Pl.c(i,:),'MarkerEdgeColor',Pl.c(i,:))
axis([Pl.xlim 0 2])
set(gca,'XTick',Pl.XTick,'XTickLabel',Pl.XTickLabel)
xlabel('Time (sec)','FontSize',Pl.fs)