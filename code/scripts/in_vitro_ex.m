clear, clc, %close all
path         = '~\Research\data\2P\rafa\adam\2009\20090128a';
imname       = '20090128a10';
tifname     = [path imname '.tif'];
paqname     = [path '20090128a10.paq'];
prettytif   = [path 'AVG_20090128a10_ROI.tif'];
figdir      = '~\Research\publications\peer-reviewed\fast-oopsi\figs\';
datadir     = '~\Research\publications\peer-reviewed\fast-oopsi\data\';
fname       = 'in_vitro_ex';

% set switches of things to do
LoadTif     = 1;
GetROI      = 1;
GetEphys    = 1;

%% get image data

if LoadTif == 1                                     % get whole movie
    MovInf  = imfinfo(tifname);                     % get number of frames
    Im.T  = numel(MovInf);                  % only alternate frames have functional data
    Im.h  = MovInf(1).Height;
    Im.w  = MovInf(1).Width;
    Im.Np = Im.w*Im.h;
    Im.fname = fname;

    Im.DataMat = zeros(Im.w*Im.h,Im.T);% initialize mat to store movie
    for j=1:Im.T
        X = imread(tifname,j);
        Im.DataMat(:,j)=X(:);
    end
    Im.MeanFrame=mean(Im.DataMat,2);
    save([datadir fname],'Im')
else
    load([datadir fname])
end

%% select roi

if GetROI == 1
    figure(1); clf,
    imagesc(reshape(Im.MeanFrame,Im.h,Im.w)')
    title('select roi radius, double click when complete')
    set(gca,'BusyAction','queu','Interruptible','on');
    ellipse0=imellipse;
    wait(ellipse0);
    vertices0=getVertices(ellipse0);
    xdat=vertices0(:,1);
    ydat=vertices0(:,2);
    Im.x0 = min(xdat) + .5*(max(xdat)-min(xdat));
    Im.y0 = min(ydat) + .5*(max(ydat)-min(ydat));
    a = max(xdat)-Im.x0;
    b = max(ydat)-Im.y0;
    Im.radius0=mean([a b]);

    [pixmatx pixmaty] = meshgrid(1:Im.h,1:Im.w);
    Im.roi = (((pixmatx-Im.x0).^2 + (pixmaty-Im.y0).^2 )<= Im.radius0^2);
    Im.roi_edge = edge(uint8(Im.roi));
    save(fname,'Im')
else
    load([datadir fname])
end

%% get spike data

if GetEphys ==1
    Ep       = paqread(paqname, 'info');
    Ep       = Ep.ObjInfo;
    channels= 1:length(Ep.Channel);
    Data    = paqread(paqname,'Channels',channels);
    for i=channels
        Ep.Channel(i).Data = Data(:,i);
    end

    j=1;
    for i=channels
        if Ep.Channel(i).ChannelName(1)=='V'
            Ep.spt{j}=GetSpikeTimes(Ep.Channel(i).Data,0.5);
            j=j+1;
        elseif strcmp(Ep.Channel(i).ChannelName,'CameraSync')
            CameraSync = Data(:,i);
        end
    end
    Ep.Nc = j-1;
    dC                  = diff(CameraSync);
    frame_onset_times   = find(dC>2);
    fake_frames         = length(frame_onset_times)-Im.T;
    frame_onset_times(1:fake_frames)=[];
    Im.dt   = median(diff(frame_onset_times))/Ep.SampleRate;
    
    Ep.n=zeros(Ep.Nc,Ep.SamplesAcquired);
    for i=1:Ep.Nc
        Ep.n(i,Ep.spt{i})=1;
    end

    Ep.n=zeros(Im.T,Ep.Nc);
    for i=1:Ep.Nc
        Ep.n(:,i) = SubSampleSpikeTrain(frame_onset_times,Ep.spt{i});
    end

    save([datadir fname],'Ep','-append')
else
    load([datadir fname])
end

%% plot ROI

Pl.xlims=[1 Im.T];
Pl.nticks=4;
Pl = PlotParams(Pl);

roi2=Im.roi';
roi_edge2=Im.roi_edge';

ROI_im      = Im.MeanFrame+max(Im.MeanFrame)*roi_edge2(:)/2;
weighted_ROI= Im.MeanFrame.*roi2(:);

fig=figure(2); clf,
height=0.4;
width=Im.h/Im.w*height;
left = (1-width)/2;
bottom=(.5-height)/2+.5;
subplot('position',[left bottom width height]);
Im.MeanPretty = imread(prettytif);
imagesc(Im.MeanPretty)
colormap('gray')
title('mean frame','FontSize',Pl.fs)

subplot(2,2,3); hold all
F=weighted_ROI'*Im.DataMat(:,1:end)/sum(weighted_ROI(:)); F=z1(F(2:end))';
plot(F+1,'k','LineWidth',Pl.lw);
x=fmincon(@(x) sum((filter(1,[1 -x(1)],x(2)*Ep.n(2:end,1))+x(3) - F).^2),[.99 1 .2],[],[],[],[],[0 0 -1], [1 1 1]);
C = filter(1,[1 -x(1)],x(2)*Ep.n(:,1))+x(3);               % calcium concentration
plot(C+1,'color',0.75*[1 1 1],'LineWidth',Pl.lw);
bar(Ep.n(:,1))
axis('tight')
xlabel('time (sec)','FontSize',Pl.fs)
set(gca,'XTick',[200:200:2000],'XTickLabel',[3 6 9 12 15],'YTickLabel',[])
ylab=ylabel([{'F'}; {''}; {''}; {''}; {''}; {''}; {'n'}],'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

subplot(2,2,4); hold all
resid=F-C(2:end);
[n,xout] = hist(resid,20);
normn=n/sum(n);
plot(xout,normn,'--k','LineWidth',Pl.lw)
ylabel('probability','FontSize',Pl.fs)
xlabel('residual error','FontSize',Pl.fs)

[muhat,sigmahat] = normfit(resid);
gauss=1/sqrt(2*pi*sigmahat^2)*exp(- (linspace(min(xout),max(xout),length(xout)) - muhat).^2/(2*sigmahat^2));
plot(xout,gauss/sum(gauss),'-k','LineWidth',Pl.lw)
axis([-.4 .4 0 max(normn)])

%if Vsave==1 % print fig
    wh=[7 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir fname];
    print('-depsc',figname)
    print('-dpdf',figname)
    saveas(fig,figname)
%end