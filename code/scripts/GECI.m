% this script simulates data with a poisson observation model, and then
% infers spikes using both the gaussian and poisson assumptions. note that
% i comment out the line that normalizes F to be betweeen eps and 1 in
% fast_oopsi because i am feeding in the parameters for this code (ie, we
% have not written code to estimate {a,b} for the poisson case).

clear, clc,
fname = 'GECI';

% set simulation metadata
T               = 10000;     % # of time steps
V.dt            = 0.001;    % time step size
V.fast_iter_max = 2;        % whether to plot with each iteration
V.fast_plot     = 1;
V.save          = 1;        % whether to save results

% initialize params
P.a     = 1/30;                % scale
P.b     = 1;           % bias
tau     = [0.5 0.03];                % decay time constant for each cell
P.gam   = 1-V.dt./tau;       % set gam
P.lam   = 1;              % rate
P.sig   = 0.15;

% simulate data
V.n = poissrnd(P.lam*V.dt*ones(T,1));   % simulate spike train
% V.n(1) = 1;
C1  = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
C2  = filter(1,[1 -P.gam(2)],-V.n);         % calcium concentration
V.C = C1+C2;
s=1;
for t=1:T
        if mod(t,30)==0;
            V.F(s) = P.a*sum(V.C(t-29:t))+P.b+P.sig*randn;
            nsub(s) = sum(V.n(t-29:t));
            s=s+1;
        end
end
figure(1), clf, hold off; plot(z1(V.F)+1); hold all; bar(V.n),

if V.save, save(['../../data/' fname '.mat'],'V','P'); end

%% infer spike trains and parameters
% V.F=z1(V.F);
V.dt=30*V.dt;
P.gam = 1 - V.dt/mean(tau);
I{1}.n = fast_oopsi(V.F,V);

if V.save, save(['../../data/' fname '.mat']); end

%% plot results
T=length(V.F);


fig=figure(1); clf, hold on
nrows   = 2;
ncols   = 1;

Pl.g    = [0.75 0.75 0.75];     % gray color
Pl.fs   = 12;                   % font size
Pl.ms   = 5;                    % marker size
Pl.lw   = 1;                    % line width
Pl.n    = nsub; Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.shift= .1;
Pl.xlim = [2 T-3];
Pl.xlims= Pl.xlim(1):Pl.xlim(2);
Pl.XTick= Pl.xlim(1):round(mean(T)/5):Pl.xlim(2);
Pl.XTickLabel = round((Pl.XTick-min(Pl.XTick))*V.dt*100)/100;

% plot F and n for neurons
subplot(nrows,ncols,1), hold on
plot((V.F),'Color','k','LineWidth',Pl.lw)
axis('tight')
% axis([Pl.xlim(1) Pl.xlim(2) 0 max(V.F)])
set(gca,'YTick',[0:10:max(V.F)]);%,'YTickLabel',[])
ylab=ylabel([{'fluorescence'}; ],'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% title('Fluorescence Projection','FontSize',Pl.fs)
set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',[]); %(Pl.XTick-min(Pl.XTick))*V.dt)
% set(gca,'YTick',[])
% xlabel('Time (sec)','FontSize',Pl.fs)

% plot inferred spike train
for q=1%:length(I)
    subplot(nrows,ncols,1+q)
    hold on
    stem(Pl.n+Pl.shift,'LineStyle','none','Marker','v','MarkerEdgeColor',Pl.g,'MarkerFaceColor',Pl.g,'MarkerSize',Pl.ms)
    bar(I{q}.n(Pl.xlims)/max(I{q}.n(Pl.xlims)),'EdgeColor','k','FaceColor','k')
    axis([Pl.xlim(1) Pl.xlim(2) 0 1+Pl.shift])
    set(gca,'YTick',[0 1],'YTickLabel',[])
    if q==1
        ylab=ylabel([{'fast'}; {'filter'}],'FontSize',Pl.fs);
        set(gca,'XTick',Pl.XTick,'XTickLabel',[])
    else
        ylab=ylabel([{'Poisson'}; {'observations'}],'FontSize',Pl.fs);
        set(gca,'XTick',Pl.XTick,'XTickLabel',(Pl.XTick-min(Pl.XTick))*V.dt)
        xlabel('Time (sec)','FontSize',Pl.fs)
    end
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
end

% print fig
if V.save
    wh=[7 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    DirName='../../figs/';
    print('-depsc',[DirName fname])
    print('-dpdf',[DirName fname])
    saveas(fig,[DirName fname])
end