% this script runs fast_oopsi and wiener_oopsi for a range of firing rates,
% to facilitate computing various statistics
clear, clc,
fname = 'stats';
cd ~/Research/oopsi/fast-oopsi/code/scripts/

% set simulation metadata
V.T             = 1000;    % # of time steps
V.dt            = 1/30;    % time step size
V.fast_iter_max = 1;        % whether to plot with each iteration
V.fast_plot     = 0;
V.save          = 0;        % whether to save results

%% vary sigma
% initialize params
P.a     = 1;                % scale
P.b     = 0;                % bias
tau     = 0.5;              % decay time constant for each cell
P.gam   = 1-V.dt./tau;      % set gam
P.lam   = 3;                % rate
sigs = 0.1:0.05:0.6;
auc.f = zeros(10,length(sigs));
auc.w = auc.f;

for j=1:length(sigs)
    P.sig=sigs(j);
    for i=1:10 % trials
        V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
        V.n(V.n>1)=1;
        V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
        V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
        fast    = fast_oopsi(V.F,V,P); 
        fast    = fast/max(fast);
        wiener  = wiener_oopsi(V.F,V.dt,P); 
        wiener  = wiener/max(abs(wiener));
        rocf    = roc3([fast V.n]);
        rocw    = roc3([wiener V.n]);
        auc.f(i,j) = rocf.AUC;
        auc.w(i,j) = rocw.AUC;
    end
end
if V.save, save(['../../data/' fname '.mat'],'V','P','auc','rocf','rocw','sigs'); end

%% get smooth roc
P.sig   = 0.35;
V.T     = 10000;
V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
V.n(V.n>1)=1;
V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
fast    = fast_oopsi(V.F,V,P);
fast    = fast/max(fast);
wiener  = wiener_oopsi(V.F,V.dt,P);
wiener  = wiener/max(abs(wiener));
roc.f    = roc3([fast V.n]);
roc.w    = roc3([wiener V.n]);


%% vary lambda
% initialize params
P.a     = 1;                % scale
P.b     = 0;                % bias
tau     = 0.5;              % decay time constant for each cell
P.gam   = 1-V.dt./tau;      % set gam
lams    = [1 3 5 10 30 50 100];                % rate
P.sig   = 0.2;
mse.f   = zeros(10,length(lams));
mse.w   = mse.f;
V.T     = 1000;
for j=1:length(lams)
    P.lam=lams(j);
    for i=1:10 % trials
        V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
        V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
        V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
        fast    = fast_oopsi(V.F,V,P);
        fast    = fast(2:end-1);
        fast    = fast/max(fast)*max(V.n);
        wiener  = wiener_oopsi(V.F,V.dt,P);
        wiener  = wiener(2:end-1);
        wiener(wiener<0) = 0;
        wiener  = wiener/max(abs(wiener))*max(V.n);
        mse.f(i,j) = sum((fast-V.n(2:end-1)).^2);
        mse.w(i,j) = sum((wiener-V.n(2:end-1)).^2);
    end
%     figure(2), clf, plot(z1(V.F)+max(V.n)), hold on, bar(V.n);  pause
end
if V.save, save(['../../data/' fname '.mat'],'roc','mse',--append); end



%% plot results
% fname = 'stats';
% load(['../../data/' fname '.mat'])
% sigs = 0.1:0.05:0.6; V.T = 10000; V.dt = 1/30; P.gam = 1-V.dt./0.5; P.a=1;P.b=0;
fig=figure(1); clf, hold on
nrows   = 2;
ncols   = 2;
fs      = 14;
xticks  = [0:90:V.T];

% plot example traces
subplot(nrows,ncols,1)
P.lam   = 3; 
P.sig   = sigs(1);
V.T     = 180;    % # of time steps
V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
V.n(V.n>1)=1;
V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
plot(z1(V.F)+2,'k')
hold all
P.sig   = sigs(end);
V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
plot(z1(V.F)+1,'k')
bar(V.n)
axis('tight')
set(gca,'XTick',xticks,'XTickLabel',xticks*V.dt,'FontSize',fs)
set(gca,'YTickLabel',[])
xlabel('time (sec)','FontSize',fs)
ylab=ylabel([{'high SNR'}; {''}; {''}; {'low SNR'}; {''}; {''}; {'spikes'}]);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',fs)

% plot ROC
subplot(nrows,ncols,3), hold all
plot(rocf.xr,rocf.yr,'-k')
plot(rocw.xr,rocw.yr,'-.','Color','k')
axis([0 0.3 0 1])
box on
set(gca,'FontSize',fs)
ylabel([{'true positive rate'}],'FontSize',fs)
xlabel('false positive rate','FontSize',fs)
legend('fast','wiener','Location','SouthEast')

% plot auc
subplot(nrows,ncols,4), cla
errorbar(mean(auc.f),std(auc.f),'-k')
hold all
errorbar(mean(auc.w),std(auc.w),'-.k')
axis('tight')
set(gca,'XTick',1:2:length(sigs),'XTickLabel',sigs(1:2:end),'FontSize',fs)
ylab=ylabel('AUC','FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlabel('standard deviation','FontSize',fs)

% plot mse
subplot(nrows,ncols,2), cla
errorbar(mean(mse.f),std(mse.f),'-k')
hold all
errorbar(mean(mse.w),std(mse.w),'-.k')
axis('tight')
set(gca,'XTick',1:2:length(lams),'XTickLabel',lams(1:2:end))
set(gca,'YTick',[1 10 100 1000 10000])
ylab=ylabel('MSE','FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlabel('expected firing rate','FontSize',fs)
set(gca,'YScale','log','FontSize',fs)
xlim([0.9 7.1])

% print fig
if V.save
    wh=[8 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    DirName='../../figs/';
    print('-depsc',[DirName fname])
    print('-dpdf',[DirName fname])
    saveas(fig,[DirName fname])
end