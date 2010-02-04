% this script runs fast_oopsi and wiener_oopsi for a range of firing rates,
% to facilitate computing various statistics
clear all; close all; clc;
fname = 'stats_data';

% change to script directory
cd ~/Research/oopsi/fast-oopsi/code/scripts/
% cd   /Work/fast-oopsi/code/scripts/

% load datasets
%dataset = load('C:\Users\Tim\Desktop\local data\various\vogelstein\Imaging-SNR-Data.mat');
% dataset  = load('/Work/local-oopsi/Imaging-SNR-Data.mat');
dataset = load('~/Research/oopsi/meta-oopsi/data/rafa/adam/2008/Imaging-SNR-Data.mat');

% this dataset has a really crazy baseline and so it isn't fair to use
dataset = rmfield(dataset,'ExpDat20080801b4');

% get names of datasets
names = fieldnames(dataset);

% set various parameters
V.fast_iter_max = 5;        % whether to plot with each iteration
V.fast_plot     = 0;        % whether to generate foopsi plots
V.save          = 0;        % whether to save results
P.a       = 1;              % scale
P.b       = 0;              % bias
dPlot     = 1;              % generate debug plot
ephysRate = 10000;          % electrophysiology sampling rate
alpha     = 0.6;            % spike detection threshold
jitters   = [1:10];
offsets   = -max(jitters):1:max(jitters);
auc.f     = zeros(length(jitters)+length(offsets),length(names));
auc.w     = auc.f;

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
    V.n = get_spike_times(vr,alpha);
    % plot spike times from ephys data
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
    % plot spike times
    if dPlot
        figure(666); ax(1) = subplot(4,1,2);
        cla; plot(V.F./max(V.F),'k'); hold on; plot(find(V.n > 0),1,'r.');
    end
    % run all oopsies
    [fast Pb]   = fast_oopsi(V.F,V);
    fast        = fast/max(fast);
    FAST{j}     = fast;
    % set up parameters for the wiener filter
    PP = P;
    if 0
        % set gam
        PP.gam = 1-V.dt./tau;
        % set lam for wiener (expected time jump per bin)
        PP.lam   = sum(fast)/(V.T*V.dt);
        % set sig for wiener (std of noise)
        PP.sig   = PP.lam;
    else
        % use foopsi estimates
        %PP = Pb;
        PP.gam = 1-V.dt/1;
        PP.sig = mad(V.F,1)*1.4826;
        PP.lam = PP.sig;
    end
    wiener  = wiener_oopsi(V.F,V.dt,PP);
    wiener  = wiener/max(abs(wiener));
    WIENER{j} = wiener;
    % plot inference output
    if dPlot
        figure(666); ax(2) = subplot(4,1,3); title(num2str(j))
        cla; bar(wiener/max(wiener),'k'); hold on;
        plot(find(V.n > 0),1,'r.');
        ax(3) = subplot(4,1,4);
        cla; bar(fast/max(fast),'k'); hold on;
        plot(find(V.n > 0),1,'r.');
        linkaxes(ax,'x');
    end
    % save roc curves
    aa = .05;
    % increase acceptable window size from 1:jitter
    for jitter = 1:length(jitters)
        rocf   = roc3_gamma([fast V.n],aa,jitters(jitter));
        rocw   = roc3_gamma([wiener V.n],aa,jitters(jitter));
        rocr   = roc3_gamma([randn(size(V.n)) V.n],aa,jitters(jitter));
        auc.f(jitter,j) = rocf.AUC;
        auc.w(jitter,j) = rocw.AUC;
        auc.r(jitter,j) = rocr.AUC;
    end
    % offset spike times by all possile offsets
    for offset = 1:length(offsets)
        rocf   = roc3_gamma([fast V.n],aa,offsets(offset),false);
        rocw   = roc3_gamma([wiener V.n],aa,offsets(offset),false);
        rocr   = roc3_gamma([randn(size(V.n)) V.n],aa,offsets(offset),false);
        auc.f(offset+length(jitters),j) = rocf.AUC;
        auc.w(offset+length(jitters),j) = rocw.AUC;
        auc.r(offset+length(jitters),j) = rocr.AUC;        
    end
    mse.f(j)=sum((fast-V.n).^2);
    temp=wiener; temp(temp<0)=0;
    mse.w(j)=sum((temp-V.n).^2);
    mse.r(j)=sum((randn(size(V.n))-V.n).^2);
end

% plot wiener auc values vs. foopsi auc values
figure; plot(rocf.xr,rocf.yr,'k.'); hold on; plot(rocw.xr,rocw.yr,'r.');
colors = jet(length(jitters));
figure; subplot(1,2,1); hold on;
for j=1:length(jitters)
    plot(auc.f(1:length(jitters),j),auc.w(1:length(jitters),j),'-.','Color',colors(j,:)); xlabel('fast'); ylabel('wiener');
end
plot([0 1],[0 1],'k'); axis square;
colors = jet(length(offsets));
subplot(1,2,2); hold on;
for j=1:length(names)
    plot(auc.f(length(jitters)+1:end,j),auc.w(length(jitters)+1:end,j),'-.','Color',colors(j,:)); xlabel('fast'); ylabel('wiener');
end

if V.save
    save(['../../data/' fname '.mat'],'V','P','auc','jitters','offsets','FAST','WIENER');
end

%% plot auc box plots
plot([0 1],[0 1],'k'); axis square;
fs=18;
figure(777); clf; subplot(2,1,1); hold on;
boxplot([auc.f(21,:)' auc.f(1:length(jitters),:)'],'plotstyle','compact','color','b','orientation','vertical');
boxplot([auc.w(21,:)' auc.w(1:length(jitters),:)'],'plotstyle','compact','color','r','orientation','vertical');
% set(gca,'YTick',0.6:0.2:1)
axis([-0.5 10.5 0.6 1.01])
xlabel('spike window size (bins)','FontSize',fs);
ylab=ylabel('AUC','FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',fs)
set(gca,'XTickLabel',0:10,'FontSize',fs);
subplot(2,1,2); hold on;
boxplot(auc.f(length(jitters)+1:end,:)','plotstyle','compact','color','b','orientation','vertical');
boxplot(auc.w(length(jitters)+1:end,:)','plotstyle','compact','color','r','orientation','vertical');
xtick=offsets(1:end);
set(gca,'XTick',1:2:21,'XTickLabel',offsets(1:2:end),'FontSize',fs)
% set(gca,'YTick',0.6:0.2:1)
xlabel('spike train jitter (bins)','FontSize',fs);
ylab=ylabel('AUC','FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',fs)

% save things if desired
if V.save
    figname=[fname num2str(1)];
    wh=[7 5];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print(gcf,'-dpdf',['../../figs/' figname '.pdf']);
    print(gcf,'-depsc2',['../../figs/' figname '.eps']);
end

%% plot auc errorbars
plot([0 1],[0 1],'k'); axis square;
fs=18;
lw=2;
gray=0.7*[1 1 1];
figure(888); clf;
subplot(2,1,1); hold all;

fmean=mean([auc.f(21,:)' auc.f(1:length(jitters),:)']);
fstd=std([auc.f(21,:)' auc.f(1:length(jitters),:)']);
errorbars(fmean,fstd,'-','color','k','linewidth',lw);

wmean=mean([auc.w(21,:)' auc.w(1:length(jitters),:)']);
wstd=std([auc.w(21,:)' auc.w(1:length(jitters),:)']);
errorbars(wmean,wstd,'-.','color',gray,'linewidth',lw);

rmean=mean([auc.r(21,:)' auc.r(1:length(jitters),:)']);
rstd=std([auc.r(21,:)' auc.r(1:length(jitters),:)']);
errorbars(wmean,wstd,'-.','color','b','linewidth',lw);

xlabel('spike window size (bins)','FontSize',fs);
set(gca,'XLim',[0.5 11.5])
ylab=ylabel('AUC','FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',fs)
set(gca,'XTick',1:11,'XTickLabel',0:10,'FontSize',fs);


subplot(2,1,2); hold all;
fmean=mean(auc.f(length(jitters)+1:end,:)');
fstd=std(auc.f(length(jitters)+1:end,:)');
errorbars(fmean,fstd,'-','color','k','linewidth',lw);
wmean=mean(auc.w(length(jitters)+1:end,:)');
wstd=std(auc.w(length(jitters)+1:end,:)');
errorbars(wmean,wstd,'-.','color',gray,'linewidth',lw);
xtick=offsets(1:end);
set(gca,'XLim',[0 22])
set(gca,'XTick',1:2:21,'XTickLabel',offsets(1:2:end),'FontSize',fs)
% set(gca,'YTick',0.6:0.2:1)
xlabel('spike train jitter (bins)','FontSize',fs);
ylab=ylabel('AUC','FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',fs)


% save things if desired
if V.save
    figname=[fname num2str(2)];
    wh=[7 5];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print(gcf,'-dpdf',['../../figs/' figname '.pdf']);
    print(gcf,'-depsc2',['../../figs/' figname '.eps']);
end

%% histograms
FastTotal=[];
WienerTotal=[];
for j=1:length(names)
    FastTotal=[FastTotal; FAST{j}];
    WienerTotal=[WienerTotal; WIENER{j}];
end

FastTotal(FastTotal<1e-5)=[];
WienerTotal(WienerTotal<1e-5)=[];
figure(1), clf,
nbins=50;
subplot(211),
[n,fout] = hist(FastTotal,nbins);
[fv fq]=getElbows(fout);
hist(FastTotal,nbins)
subplot(212),
[n,wout] = hist(WienerTotal,nbins);
[wv wq]=getElbows(wout);
hist(WienerTotal,nbins)

fq(end-1:end)=[];
wq(end-1:end)=[];

%% plot mse

for j=1:length(names)
    cc = dataset.(char(names(j)));
    V.F = cc.Fluorescence;
    vr = cc.chanDev1_ai0_VoltageCh1;
    V.n = get_spike_times(vr,0.6);
    for ii = 1:length(V.n)
        [value V.n(ii)] = min(abs(cc.FluorescenceTime*ephysRate - V.n(ii)));
    end
    times = zeros(length(V.F),1);
    times(V.n) = 1; V.n = times;

    temp=FAST{j};
    temp=temp/max(temp)*max(V.n);
    mse.f(j)=sum((temp-V.n).^2);
    temp=WIENER{j};
    temp(temp<0)=0;
    temp=temp/max(temp)*max(V.n);
    mse.w(j)=sum((temp-V.n).^2);
    
    mse.c(j)=sum((0*ones(size(V.n))-V.n).^2);
end

fig=figure(999); clf, hold all
H=boxplot([mse.f; mse.w; mse.c]','notch','on');
set(H([1 2 3 4 5 6],:),'LineWidth',2,'Color','k')
set(H([1 2],:),'LineStyle','-')
set(H(7,2),'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',2)
set(gca,'XTickLabel',[{'fast'}; {'wiener'}; {'zeros'}],'FontSize',fs)
xlabel([])
ylabel('mean-squared error')
% set(gca,'YScale','log')
% set(gca,'YLim',[0 100])

% save things if desired
if V.save
    figname=['boxplots_data'];
    wh=[7 5];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print(gcf,'-dpdf',['../../figs/' figname '.pdf']);
    print(gcf,'-depsc2',['../../figs/' figname '.eps']);
end

