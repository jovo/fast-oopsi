% this script runs fast_oopsi and wiener_oopsi for a range of firing rates,
% to facilitate computing various statistics
clear, clc,
fname = 'stats_data';

% change to script directory
%cd ~/Research/oopsi/fast-oopsi/code/scripts/
cd C:\Users\Tim\oopsi\fast-oopsi\fast-oopsi\code\scripts\

% load datasets
dataset = load('C:\Users\Tim\Desktop\local data\various\vogelstein\Imaging-SNR-Data.mat');

% get names of datasets
names = fieldnames(dataset);

% set various parameters
V.fast_iter_max = 5;        % whether to plot with each iteration
V.fast_plot     = 0;        % whether to generate foopsi plots
V.save          = 0;        % whether to save results
P.a       = 1;              % scale
P.b       = 0;              % bias
tau       = 0.5;            % decay time constant for each cell
auc.f     = zeros(1,length(names));
auc.w     = auc.f;
dPlot     = 1;              % generate debug plot
ephysRate = 10000;          % electrophysiology sampling rate

% estimate all parameters
est_sig = 1;  
est_lam = 1;
est_gam = 1;
est_b   = 1;
est_a   = 1;

%% iterate over real datasets
for j=1:length(names)
    % get a dataset struct
    cc = dataset.(char(names(j)));
    % compute frame rate
    V.dt =mean(diff(cc.FluorescenceTime));
    % get the fluorescence trace
    V.F = cc.Fluorescence;
    % normalize the trace
    V.F=V.F-min(V.F); V.F=V.F/std(V.F);
    V.F=V.F/max(V.F); V.F=V.F+realmin;
    % set parameters based on real data
    V.T = length(V.F);
    % get the intracellularly measured spike train
    vr = cc.chanDev1_ai0_VoltageCh1;
    V.n = get_spike_times(vr);
    % make sure we're not stupid
    if dPlot
        figure(666); subplot(4,1,1);
        cla; plot(vr./max(vr),'k');
        hold on; plot(V.n,1,'r.');
    end
    % convert spike times from ephys samples to movie frames
    for ii = 1:length(V.n)
        [value V.n(ii)] = min(abs(cc.FluorescenceTime*ephysRate - V.n(ii)));
    end
    times = zeros(length(V.F),1);
    times(V.n) = 1; V.n = times;
    % plot spike times for fun
    if dPlot
        figure(666); subplot(4,1,2);
        cla; plot(V.F./max(V.F),'k'); hold on; plot(V.n,'r.');
    end
    % run all oopsies
    [fast Pb]    = fast_oopsi(V.F,V,P);
    fast    = fast/max(fast);
    % use below estimates (good) or foopsi estimates (bad)
    PP = P;
    if 1
        % set gam
        PP.gam = 1-V.dt./tau;
        % set lam for wiener
        PP.lam   = sum(fast)/(V.T*V.dt);
        % set sig for wiener
        PP.sig   = PP.lam;
    else
        % use foopsi estimates
        PP = Pb;
    end
    wiener  = wiener_oopsi(V.F,V.dt,PP);
    wiener  = wiener/max(abs(wiener));
    % plot inference output
    if dPlot
        figure(666); subplot(4,1,3);
        cla; bar(wiener/max(wiener),'k'); hold on;
        plot(V.n,'r.');
        subplot(4,1,4);
        cla; bar(fast/max(fast),'k'); hold on; 
        plot(V.n,'r.');
        keyboard;
    end
    % save roc curves
    rocf    = roc3([fast V.n]);
    rocw    = roc3([wiener V.n]);
    auc.f(j) = rocf.AUC;
    auc.w(j) = rocw.AUC;
    % save things if desired
    if V.save, save(['../../data/' fname '.mat'],'V','P','auc'); end
end




% %% get smooth roc
% P.sig   = 0.35;
% V.T     = 10000;
% V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
% V.n(V.n>1)=1;
% V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
% V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
% fast    = fast_oopsi(V.F,V,P);
% fast    = fast/max(fast);
% wiener  = wiener_oopsi(V.F,V.dt,P);
% wiener  = wiener/max(abs(wiener));
% roc.f    = roc3([fast V.n]);
% roc.w    = roc3([wiener V.n]);
% 
% 
% %% vary lambda
% % initialize params
% P.a     = 1;                % scale
% P.b     = 0;                % bias
% tau     = 0.5;              % decay time constant for each cell
% P.gam   = 1-V.dt./tau;      % set gam
% lams    = [1 3 5 10 30 50 100];                % rate
% P.sig   = 0.2;
% mse.f   = zeros(10,length(lams));
% mse.w   = mse.f;
% V.T     = 1000;
% for j=1:length(lams)
%     P.lam=lams(j);
%     for i=1:10 % trials
%         V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
%         V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
%         V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
%         fast    = fast_oopsi(V.F,V,P);
%         fast    = fast(2:end-1);
%         fast    = fast/max(fast)*max(V.n);
%         wiener  = wiener_oopsi(V.F,V.dt,P);
%         wiener  = wiener(2:end-1);
%         wiener(wiener<0) = 0;
%         wiener  = wiener/max(abs(wiener))*max(V.n);
%         mse.f(i,j) = sum((fast-V.n(2:end-1)).^2);
%         mse.w(i,j) = sum((wiener-V.n(2:end-1)).^2);
%     end
% %     figure(2), clf, plot(z1(V.F)+max(V.n)), hold on, bar(V.n);  pause
% end
% if V.save, save(['../../data/' fname '.mat'],'auc','roc','mse'); end



%% plot results
% fig=figure(1); clf, hold on
% nrows   = 2;
% ncols   = 2;
% fs      = 14;
% xticks  = [0:90:V.T];
% 
% % plot example traces
% subplot(nrows,ncols,1)
% P.lam   = 3; 
% P.sig   = sigs(1);
% V.T     = 180;    % # of time steps
% V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
% V.n(V.n>1)=1;
% V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
% V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
% plot((V.F./max(V.F))+2,'k')
% hold all
% P.sig   = sigs(end);
% V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
% plot((V.F./max(V.F))+1,'k')
% bar(V.n)
% axis('tight')
% set(gca,'XTick',xticks,'XTickLabel',xticks*V.dt,'FontSize',fs)
% set(gca,'YTickLabel',[])
% xlabel('time (sec)','FontSize',fs)
% ylab=ylabel([{'high SNR'}; {''}; {''}; {'low SNR'}; {''}; {''}; {'spikes'}]);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',fs)
% 
% % plot ROC
% subplot(nrows,ncols,3), hold all
% plot(rocf.xr,rocf.yr,'-k')
% plot(rocw.xr,rocw.yr,'-.','Color','k')
% axis([0 0.3 0 1])
% box on
% set(gca,'FontSize',fs)
% ylabel([{'true positive rate'}],'FontSize',fs)
% xlabel('false positive rate','FontSize',fs)
% legend('fast','wiener','Location','SouthEast')
% 
% % plot auc
% subplot(nrows,ncols,4), cla
% errorbar(mean(auc.f),var(auc.f),'-k')
% hold all
% errorbar(mean(auc.w),var(auc.w),'-.k')
% axis('tight')
% set(gca,'XTick',1:2:length(sigs),'XTickLabel',sigs(1:2:end),'FontSize',fs)
% ylab=ylabel('AUC','FontSize',fs);
% % set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlabel('standard deviation','FontSize',fs)
% 
% % plot mse
% subplot(nrows,ncols,2), cla
% errorbar(mean(mse.f),std(mse.f),'-k')
% hold all
% errorbar(mean(mse.w),std(mse.w),'-.k')
% axis('tight')
% set(gca,'XTick',1:2:length(lams),'XTickLabel',lams(1:2:end))
% set(gca,'YTick',[1 10 100 1000 10000])
% ylab=ylabel('MSE','FontSize',fs);
% % set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlabel('expected firing rate','FontSize',fs)
% set(gca,'YScale','log','FontSize',fs)
% xlim([0.9 7.1])
% 
% % print fig
% if V.save
%     wh=[8 6];   %width and height
%     set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
%     DirName='../../figs/';
%     print('-depsc',[DirName fname])
%     print('-dpdf',[DirName fname])
%     saveas(fig,[DirName fname])
% end