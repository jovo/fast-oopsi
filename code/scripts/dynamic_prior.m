% this script is a wrapper that loads some data (or simulates it) and then
% calls the most recent FOOPSI function to infer spike trains and
% parameters

clear, clc,
fname = 'dynamic_prior';

%% set parameters

% set meta parameters
Meta.T      = 600;                              % # of time steps
Meta.dt     = 1/60;                             % time step size
Meta.Np     = 1;                                % # of pixels in each image
Meta.Nc     = 1;                                % # cells per ROI
Meta.Poiss  = 0;                                % whether observations are poisson or gaussian
Meta.MaxIter= 25;                               % # iterations of EM to estimate params

% initialize other parameters
P.a     = 1;
P.b     = 0.01*mean(P.a);                       % baseline is a scaled down version of the sum of the spatial filters
P.sig   = max(P.a)*14;                         % stan dev of noise (indep for each pixel)
C_0     = 0;                                    % initial calcium
tau     = rand(Meta.Nc,1)/2+.05;                % decay time constant for each cell
P.gam   = 1-Meta.dt./tau(1:Meta.Nc);            % set gam
P.lam   = (2*(sin(linspace(0,4*pi,Meta.T))+1)'); % rate-ish, ie, lam*dt=# spikes per second

%% simulate data

% simulate
n=zeros(Meta.T,Meta.Nc);                        % pre-allocate memory for spike train
C=n;                                            % pre-allocate memory for calcium
for i=1:Meta.Nc
    n(1,i)      = C_0;                              
    n(2:end,i)  = poissrnd(P.lam(2:end,i)*Meta.dt); % simulate spike train
    n(n>1)      = 1;                                % make only 1 spike per bin
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));   % calcium concentration
end
F = P.a*(C+repmat(P.b,Meta.T,1))'+P.sig*rand(Meta.Np,Meta.T);

% estimating parameters (these are only used if MaxIter>1)
Est.sig     = 0;
Est.lam     = 0;
Est.gam     = 0;
Est.b       = 0;
Est.a       = 0;

Est.Thresh  = 1;                            % whether to threshold spike train before estimating 'a' and 'b' (we always keep this on)
Est.Plot    = 0;                            % whether to plot filter with each iteration
Est.n       = n;                            % keep true spike times to also plot for comparison purposes

save(['../../data/' fname '.mat'],'F','n','P','Meta','Est')

%% infer spike trains and parameters

q=0;
exps=1:2;
for qq=exps
    disp(qq)
    if qq==1;                                % use true params
        PP=P;
        FF=F;
        Meta.MaxIter = 1;
    elseif qq==2                             % use static prior
        PP=P;
        PP.lam = mean(P.lam)+0*P.lam;
        FF=F;
        Meta.MaxIter = 1;
    end
    q=q+1;
    [I{q}.n I{q}.P] = FOOPSI_v3_05_01(FF,PP,Meta,Est);
end

save(['../../data/' fname '.mat'],'-append','I')
sound(1*sin(linspace(0,180*pi,2000)))

%% plot results

clear Pl
nrows   = 2+numel(exps);                        % set number of rows
h       = zeros(nrows,1);                       % pre-allocate memory for handles
Pl.xlims= [1 Meta.T];                           % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 4;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = Meta.Nc;

fnum = figure(1); clf,

% plot fluorescence data
i=1; h(i) = subplot(nrows,1,i);
Pl.label = [{'Fluorescence'}];
Pl.color = 'k';
plot(z1(F(Pl.xlims(1):Pl.xlims(2))),'k','LineWidth',Pl.lw);
ylab=ylabel('Fluorescence','Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','Interpreter',Pl.interp)
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims-Pl.xlims(1) 0 1])
box off

% plot inferred spike trains
for q=exps
    i=i+1; h(i) = subplot(nrows,1,i); hold on
    if q==1
        Pl.label = [{'Dynamic'}; {'Prior'}];
    else
        Pl.label = [{'Static'}; {'Prior'}];        
    end
    stem(z1(I{q}.n(Pl.xlims(1):Pl.xlims(2))),'Marker','none','LineWidth',Pl.sw,'Color','k')
    stem(Pl.n(Pl.xlims(1):Pl.xlims(2))*1.1,'Marker','v','MarkerSize',Pl.vs,...   % plot real spike train
        'LineStyle','none','MarkerFaceColor',Pl.gray,'MarkerEdgeColor',Pl.gray);
    axis([Pl.xlims min(min(n),0) 1.15])
    ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs,'Interpreter',Pl.interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:1,'YTickLabel',[])
    set(gca,'XTick',Pl.XTicks,'XTickLabel',[])
    set(gca,'XTickLabel',[])
    box off
end

% plot prior
i=i+1; h(i) = subplot(nrows,1,i);
Pl.color = Pl.colors(1,:);
plot(z1(P.lam(Pl.xlims(1):Pl.xlims(2))),'k','LineWidth',Pl.lw);
ylab=ylabel('Prior','Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','Interpreter',Pl.interp)
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims-Pl.xlims(1) 0 1.1])
box off

% set xlabel stuff
subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Meta.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
linkaxes(h,'x')

% print fig
wh=[7.5 3.5*length(exps)];   %width and height
DirName = '../../figs/';
PrintFig(wh,DirName,fname);
