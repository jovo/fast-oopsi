clear, clc, %close all
path         = '/Users/joshyv/Research/oopsi/meta-oopsi/data/rafa/adam/2009/20090128a/';
imname       = '20090128a10_cropped';
tifname     = [path imname '.tif'];
paqname     = [path '20090128a10.paq'];
fname       = 'in_vitro_ex';

% set switches of things to do
LoadTif     = 0;
GetROI      = 0;
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
    save(fname,'Im')
else
    load(fname)
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
    load(fname)
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

    save(fname,'Ep','-append')
else
    load(fname)
end

%% plot ROI
Pl = PlotParams;
roi2=Im.roi';
roi_edge2=Im.roi_edge';

ROI_im      = Im.MeanFrame+max(Im.MeanFrame)*roi_edge2(:);
weighted_ROI= Im.MeanFrame.*roi2(:);

figure(2); clf,
subplot(2,2,1:2);
imagesc(reshape(ROI_im,Im.h,Im.w)')

subplot(2,2,3); hold all
F=weighted_ROI'*Im.DataMat(:,1:end)/sum(weighted_ROI(:)); F=z1(F(2:end))';
plot(F+1,'k','LineWidth',Pl.lw);
x=fmincon(@(x) sum((filter(1,[1 -x(1)],x(2)*Ep.n(2:end,1))+x(3) - F).^2),[.99 1 .2],[],[],[],[],[0 0 -1], [1 1 1])
C = filter(1,[1 -x(1)],x(2)*Ep.n(:,1))+x(3);               % calcium concentration
plot(C+1,'r','LineWidth',Pl.lw);
bar(Ep.n(:,1))
axis('tight')

subplot(2,2,4); hold all
resid=F-C(2:end);
[n,xout] = hist(resid,20);
normn=n/sum(n);
plot(xout,normn,'k','LineWidth',Pl.lw)

[muhat,sigmahat] = normfit(resid);
gauss=1/sqrt(2*pi*sigmahat^2)*exp(- (linspace(min(xout),max(xout),length(xout)) - muhat).^2/(2*sigmahat^2));
plot(xout,gauss/sum(gauss),'Color',Pl.gray,'LineWidth',Pl.lw)
axis([-.4 .4 0 max(normn)])
