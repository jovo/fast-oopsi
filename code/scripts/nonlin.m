% this script is a wrapper that loads some data (or simulates it) and then
% calls the most recent FOOPSI function to infer spike trains and
% parameters

clear, clc,
fname = 'nonlin';

% set simulation metadata
T               = 3000;     % # of time steps
V.dt            = 1/200;    % time step size
V.fast_iter_max = 1;        % # iterations of EM to estimate params
V.fast_plot     = 1;        % whether to plot with each iteration
V.save          = 0;        % whether to save results

% initialize params
P.a     = 5;                % scale
P.b     = 0.1;              % bias
tau     = 1;                % decay time constant for each cell
P.gam   = 1-V.dt/tau;       % set gam
P.lam   = 1.0;              % rate
P.sig   = 0.01;             % standard deviation
P.k_d   = 200;              % dissociation constant

% simulate data
V.n = poissrnd(P.lam*V.dt*ones(T,1));   % simulate spike train
V.C = filter(1,[1 -P.gam],V.n);         % calcium concentration
V.S = V.C./(V.C+P.k_d);
V.F = P.a*(V.S+repmat(P.b,T,1))'+P.sig*randn(1,T);

figure(1), clf, hold off; plot(z1(V.F)+1); hold all; bar(V.n),

if V.save, save(['../../data/' fname '.mat'],'V','P'); end

%% infer spike trains and parameters

for q=1:2
    if q==1,
        V.fast_nonlin=0;
        V.fast_iter_max=10;
        %         I{q}.n=V.n;
        [I{q}.n I{q}.P I{q}.V] = fast_oopsi(V.F,V);
    else
        V.fast_nonlin=1;
        V.gauss_n=I{1}.n;
        V.fast_iter_max=1;
        [I{q}.n I{q}.P I{q}.V] = fast_oopsi(V.F,V,P);
    end
end

if V.save, save(['../../data/' fname '.mat'],'-append','I'); end

%% plot results

fig=figure(1); clf, hold on
nrows   = 1+length(I);
ncols   = 1;

Pl.g    = [0.75 0.75 0.75];     % gray color
Pl.fs   = 12;                   % font size
Pl.ms   = 5;                    % marker size
Pl.lw   = 1;                    % line width
Pl.n    = V.n; Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.shift= .1;
Pl.xlim = [2 T-3];
Pl.xlims= Pl.xlim(1):Pl.xlim(2);
Pl.XTick= Pl.xlim(1):round(mean(T)/5):Pl.xlim(2);
Pl.XTickLabel = round((Pl.XTick-min(Pl.XTick))*V.dt*100)/100;

% plot F and n for neurons
subplot(nrows,ncols,1), hold on
plot(V.F,'Color','k','LineWidth',Pl.lw)
axis([Pl.xlim(1) Pl.xlim(2) min(V.F) max(V.F)])
set(gca,'YTick',[0:10:max(V.F)]);%,'YTickLabel',[])
ylab=ylabel([{'Fluorescence'}; {'(photons/frame)'}],'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
title('Fluorescence Projection','FontSize',Pl.fs)
set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',[]); %(Pl.XTick-min(Pl.XTick))*V.dt)
% xlabel('Time (sec)','FontSize',Pl.fs)

% plot inferred spike train
for q=1:length(I)
    subplot(nrows,ncols,1+q)
    hold on
    stem(Pl.n+Pl.shift,'LineStyle','none','Marker','v','MarkerEdgeColor',Pl.g,'MarkerFaceColor',Pl.g,'MarkerSize',Pl.ms)
    bar(I{q}.n(Pl.xlims)/max(I{q}.n(Pl.xlims)),'EdgeColor','k','FaceColor','k')
    axis([Pl.xlim(1) Pl.xlim(2) 0 1+Pl.shift])
    set(gca,'YTick',[0 1],'YTickLabel',[])
    if q==1
        ylab=ylabel([{'Gaussian'}; {'Observations'}],'FontSize',Pl.fs);
        set(gca,'XTick',Pl.XTick,'XTickLabel',[])
    else
        ylab=ylabel([{'Poisson'}; {'Observations'}],'FontSize',Pl.fs);
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