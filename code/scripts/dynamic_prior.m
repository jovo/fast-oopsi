% this script is a wrapper that loads some data (or simulates it) and then
% calls the most recent FOOPSI function to infer spike trains and
% parameters

clear, clc,
fname = 'dynamic_prior';

%% set algorithmic Variables

% set meta parameters
V.T         = 1500;              % # of time steps
V.dt        = 1/60;             % time step size
V.fast_plot = 0;                % whether to plot filter with each iteration
V.fast_iter_max = 1;            % # iterations of EM to estimate params

%% initialize model Parameters
P.a     = 1;                    % scale of fluorescence data
P.b     = 0.01*mean(P.a);       % baseline is a scaled down version of the sum of the spatial filters
P.sig   = 0.5; %3*max(P.a);          % stan dev of noise (indep for each pixel)
C_0     = 0;                    % initial calcium
tau     = 1;                    % decay time constant for each cell
P.gam   = 1-V.dt./tau;          % set gam
% P.lam   = exp(2*(sin(linspace(0,10*pi,V.T))+1)')/25+1; % rate-ish, ie, lam*dt=# spikes per second
P.lam   = 0.1*ones(V.T,1); % rate-ish, ie, lam*dt=# spikes per second
spts = [1/V.dt:4/V.dt:V.T-5];
P.lam(spts)=100;
figure(2), plot(P.lam)

%% simulate data

n = poissrnd(P.lam*V.dt)';    % simulate spike train
n(n>1)=1;
C = filter(1,[1 -P.gam],n);   % convolve with exponential to get calcium
F = P.a*C+P.b+P.sig*randn(1,V.T);

save(['../../data/' fname '.mat'],'F','n','P','V')

%% infer spike trains and parameters

for q=1:2
    if q==1;                                % use true params
        PP=P;
        PP.lam=mean(P.lam);
        V.fast_iter_max = 1;
    elseif q==2                             % use static prior
        PP=P;
        V.fast_iter_max = 1;
    end
    [I{q}.n I{q}.P] = fast_oopsi(F,V,PP);
end

save(['../../data/' fname '.mat'],'-append','I')
% sound(1*sin(linspace(0,180*pi,2000)))

%% plot results

clear Pl
nrows   = 2+q;                        % set number of rows
h       = zeros(nrows,1);                       % pre-allocate memory for handles
Pl.xlims= [1 V.T];                              % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 4;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];

fnum = figure(1); clf,

% plot fluorescence data
i=1; h(i) = subplot(nrows,1,i);
Pl.label = [{'fluorescence'}];
Pl.color = 'k';
plot(z1(F(Pl.xlims(1):Pl.xlims(2))),'k','LineWidth',Pl.lw);
ylab=ylabel('fluorescence','Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','Interpreter',Pl.interp)
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims-Pl.xlims(1) 0 1])
box off

% plot inferred spike trains
for q=1:2 %exps
    i=i+1; h(i) = subplot(nrows,1,i); hold on
    if q==1
        Pl.label = [{'static'}; {'prior'}];
    else
        Pl.label = [{'dynamic'}; {'prior'}];        
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
plot((P.lam(Pl.xlims(1):Pl.xlims(2))),'k','LineWidth',Pl.lw);
ylab=ylabel('prior','Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','Interpreter',Pl.interp)
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
axis([Pl.xlims 0 max(P.lam)*1.1])
box off

% set xlabel stuff
subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*V.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
linkaxes(h,'x')

% print fig
wh=[7.5 3.5*q];   %width and height
DirName = '../../figs/';
PrintFig(wh,DirName,fname);
