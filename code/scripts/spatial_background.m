% this script generates a simulation of a movie containing a single cell
% using the following generative model:
%
% F_t = \sum_i a_i*C_{i,t} + b + sig*eps_t, eps_t ~ N(0,I)
% C_{i,t} = gam*C_{i,t-1} + n_{i,t},      n_{i,t} ~ Poisson(lam_i*dt)
%
% where ai,b,I are p-by-q matrices.
% we let b=0 and ai be the difference of gaussians (yielding a zero mean
% matrix)
%

clear, clc

% 1) generate spatial filters

% stuff required for each spatial filter
Nc      = 1;                        % # of cells in the ROI
neur_w  = 15;                       % width per neuron
width   = 15;                       % width of frame (pixels)
height  = Nc*neur_w;                % height of frame (pixels)
Npixs   = width*height;             % # pixels in ROI
x       = linspace(-5,5,height);
y       = linspace(-5,5,width);
[X,Y]   = meshgrid(x,y);
g1      = zeros(Npixs,Nc);
g2      = 0*g1;
Sigma1  = diag([1,1])*2;            % var of positive gaussian
Sigma2  = diag([1,1])*2.5;          % var of negative gaussian
mu      = [0 0];                    % means of gaussians for each cell (distributed across pixel space)
w       = [1 0.5];                  % weights of each filter


ix0 = [97:99 112:114 127:129];      % typical spatial filter
for i=1:Nc                          % spatial filter
    g1(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma1);
    g2(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma2);
end

% 2) set simulation metadata
V.T     = 400;                    % # of time steps
V.dt    = 0.005;                  % time step size
V.Np    = Npixs;                  % # of pixels in each image
V.w     = width;                  % width of frame (pixels)
V.h     = height;                 % height of frame (pixels)
V.Nc    = Nc;                     % # cells
V.plot  = 0;                      % whether to plot filter with each iteration
V.save  = 0;
V.name  = 'spatial_background';

% 3) initialize params
for i=1:V.Nc
    P.a(:,i)=g1(:,i)-g2(:,i);
end
% figure(2), imagesc(reshape(P.a,V.w,V.h)), colormap(gray), colorbar

P.b     = 0*P.a(:,1);             % baseline is zero
P.sig   = 0.01;                    % stan dev of noise (indep for each pixel)
C_0     = 0;                        % initial calcium
tau     = round(100*rand(V.Nc,1))/100+0.05;   % decay time constant for each cell
P.gam   = 1-V.dt./tau(1:V.Nc);
P.lam   = 5;                        % rate-ish, ie, lam*dt=# spikes per second

% 4) simulate data
n=zeros(V.T,V.Nc);
C=n;
for i=1:V.Nc
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(i)*V.dt*ones(V.T-1,1));    % simulate spike train
    n(n>1)      = 1;
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));           % calcium concentration
end
Z = 0*n(:,1);
F = C*P.a' + repmat(P.b',V.T,1) + P.sig*randn(V.T,Npixs);                  % fluorescence

%% 5) make movie

MakMov  = 0;
% make movie of raw data
if MakMov==1
    for i=1:V.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(F(i,:),width,height),'../../data/spatial_background_mov.tif','tif','Compression','none','WriteMode',mod)
    end
end

%% 6) infer spike train using various approaches

qs=1:2;
for q=qs
    GG=F; Tim=V;
    if q==1,
        Phat{q}=P;
        fast{q}.label='true filter';
    elseif q==2
        Phat{q}=P;
        boxcar=zeros(V.w,V.h);
        boxcar(6:10,6:10)=1;
        Phat{q}.a=boxcar(:);
        fast{q}.label='boxcar filter';
    elseif q==3
        Phat{q}=P;
        Phat{q}.a=mean(F)';
        fast{q}.label='mean filter';
    end
    display(fast{q}.label)
    starttime=cputime;
    [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GG',Tim,Phat{q});
    fast{q}.V.time = cputime-starttime;
end
if V.save==1, save(['../../data/' V.name]); end

%% plot results
% load(['../../data/' V.name])

nrows   = 3;                                 % set number of rows
ncols   = 2;
h       = zeros(nrows,1);
xlims= [5 V.T-5];                            % time steps to plot
nticks=5;                                    % number of ticks along x-axis
fs      = 14;                       % font size
ms      = 5;                        % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
Nc   = V.Nc;
XTicks=100:100:400;
    gray    = 0.4*[1 1 1];            % define gray color

fig = figure(1); clf,

ntemp=double(n);
ntemp(ntemp==0)=NaN;

for q=qs
    if q==1, i=0; else i=1; end
    i=i+1;
    subplot(nrows,ncols,i)
    imagesc(reshape(Phat{q}.a,V.w,V.h))
    colormap('gray')
    title(fast{q}.label,'FontSize',fs)
    
    i=i+2;
    subplot(nrows,ncols,i)
    plot(z1(Phat{q}.a\F'),'k','LineWidth',lw)
    hold on
    stem(ntemp,'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)

    set(gca,'XTick',XTicks,'XTickLabel',[],'FontSize',fs)
    set(gca,'YTick',[],'FontSize',fs)
    axis([0 V.T 0 1])
    if q==1, ylabel('fluorescence','FontSize',fs); end
    box off

    i=i+2;
    subplot(nrows,ncols,i)
    bar(fast{q}.n)
    hold on
    stem(ntemp,'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    axis('tight')
    set(gca,'XTick',XTicks,'XTickLabel',XTicks*V.dt,'FontSize',fs)
    set(gca,'YTick',[],'FontSize',fs)
    xlabel('time (sec)','FontSize',fs)
    if q==1, ylabel('fast filter','FontSize',fs); end
    box off
end


if V.save==1 % print fig
    wh=[7 5];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=['../../figs/' V.name];
    print('-depsc',figname)
    print('-dpdf',figname)
    saveas(fig,figname)
end