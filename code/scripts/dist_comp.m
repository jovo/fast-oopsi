%%
% N   = 1e4;
% 
% lam = 1;
% npoiss  = poissrnd(lam,N,1);
% nexp    = exprnd(lam,N,1);
% ngauss  = normrnd(lam,lam,N,1);
% 
% nbins=20;
% [nhistpoiss xout] = hist(npoiss);
% [nhistexp xout2] = hist(nexp,nbins);
% [nhistgauss xout3] = hist(ngauss,nbins);
% 
% %%
% Pl = PlotParams;
% figure(1), clf
% subplot(121)
% plot(xout,nhistpoiss,'-k','LineWidth',2), hold all
% plot(xout2, nhistexp,':k','LineWidth',2)
% plot(xout3, nhistgauss,'-.','Color',[.25 .25 .25],'LineWidth',2), hold all
% 
% 
% %%
% lam = 50;
% npoiss  = poissrnd(lam,N,1);
% nexp    = exprnd(lam,N,1);
% ngauss  = normrnd(lam,lam,N,1);
% 
% nbins=20;
% [nhistpoiss xout] = hist(npoiss);
% [nhistexp xout2] = hist(nexp,nbins);
% [nhistgauss xout3] = hist(ngauss,nbins);
% 
% %%
% subplot(122), cla
% plot(xout,nhistpoiss,'-k','LineWidth',2), hold all
% plot(xout2,nhistexp,':k','LineWidth',2)
% plot(xout3,nhistgauss,'-.','Color',[.25 .25 .25],'LineWidth',2), hold all
%%
clc, clear all, 


xmin=0;
xmax=10;
x=xmin:.01:xmax;
xint=xmin:xmax;
x2=unique([-x x]);

lam=1;
for i=xint
    xfact(i+1)=factorial(xint(i+1));
end
Npoiss=exp(-lam)*lam.^xint./xfact;

Nexp=lam*exp(-x*lam);

Ngauss=1/sqrt(2*pi*lam^2)*exp(-(x2-lam).^2/(lam^2));


figure(2), clf, 
Pl=PlotParams;
subplot(121), hold
plot(xint,Npoiss,'-k','LineWidth',2)
plot(x,Nexp,'--','Color','k','LineWidth',2)
plot(x2,Ngauss,'--','Color',.6*ones(3,1),'LineWidth',2)
plot([0 0], [0 1], '-k', 'LineWidth',0.5)
ymax=.5;
axis([-1 10 0 ymax])
title('Slow Firing Rate','FontSize',Pl.fs)
% set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',(Pl.XTick-min(Pl.XTick))*Meta.dt)
xlabel('# of spikes/frame','FontSize',Pl.fs)
ylab=ylabel('Probability','Interpreter',Pl.inter,'FontSize',Pl.fs);
% set(ylab,'Rotation',90,'HorizontalAlignment','right','verticalalignment','middle','Interpreter',Pl.interp,'FontName',Pl.font)
set(gca,'YTick',linspace(0,ymax,4),'YTickLabel',[],'TickLength',[.02 0])

%%
lam=10;
xmin=0;
xmax=30;
x=xmin:.01:xmax;
xint=xmin:xmax;
x2=unique([-x x]);

for i=xint
    xfact(i+1)=factorial(xint(i+1));
end
Npoiss=exp(-lam)*lam.^xint./xfact;

Nexp=lam*exp(-x*lam);

Ngauss=1/sqrt(2*pi*lam)*exp(-(x2-lam).^2/(lam));

subplot(122), hold
plot(x,Nexp,'--','Color','k','LineWidth',2)
plot(x2,Ngauss,'-.','Color',.6*ones(3,1),'LineWidth',2)
plot(xint,Npoiss,'-k','LineWidth',2)
plot([0 0], [0 1], '-k', 'LineWidth',0.5)
ymax=.11;
axis([0 xmax 0 ymax*1.5])
ylab=ylabel('Probability','Interpreter',Pl.inter,'FontSize',Pl.fs);
title('Fast Firing Rate','FontSize',Pl.fs)
% set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',(Pl.XTick-min(Pl.XTick))*Meta.dt)
xlabel('# of spikes/frame','FontSize',Pl.fs)
set(gca,'YTick',linspace(0,ymax,4),'YTickLabel',[],'TickLength',[.02 0])

%%
% print fig
wh=[8 3];   %width and height
DirName = '../../figs/';
PrintFig(wh,DirName,'dist_comp');
