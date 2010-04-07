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

clear,clc

% 1) generate spatial filters

% % stuff required for each spatial filter
Nc      = 1;                                % # of cells in the ROI
neur_w  = 1;                               % width per neuron
width   = 1;                               % width of frame (pixels)
height  = Nc*neur_w;                        % height of frame (pixels)
Npixs   = width*height;                     % # pixels in ROI

% 2) set simulation metadata
Sim.T       = 800;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.MaxIter = 0;                                % # iterations of EM to estimate params
Sim.Np      = Npixs;                            % # of pixels in each image
Sim.w       = width;                            % width of frame (pixels)
Sim.h       = height;                           % height of frame (pixels)
Sim.Nc      = Nc;                               % # cells
Sim.plot    = 0;                                % whether to plot filter with each iteration

lam         = [10; 500];
sigs        = [1/4 8];
% 3) initialize params
P.a     = 1;
P.b     = 0;                           % baseline is zero

P.sig   = 0.25;                                 % stan dev of noise (indep for each pixel)
C_0     = 0;                                    % initial calcium
tau     = [.1 .5]; %round(100*rand(Sim.Nc,1))/100+0.05;   % decay time constant for each cell
P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
dim     = 5;
stim    = rand(Sim.T,dim)'*2;
P.lam   = sin(linspace(0,pi,dim))';


% 3) simulate data
qs=1:10;
for q=qs
    n=zeros(Sim.T,Sim.Nc);
    C=n;
    n(1)      = C_0;
    n(2:end)  = poissrnd(P.lam'*stim(:,2:end));    % simulate spike train
    C         = filter(1,[1 -P.gam],n);               % calcium concentration
    F = C*P.a' + P.b'+P.sig*randn(Sim.T,Npixs);               % fluorescence

%% infer spikes

    D{q}.n=n; D{q}.C=C; D{q}.F=F;
    GG=D{q}.F; Tim=Sim;
    Phat{q}=P;
    Phat{q}.lam=mean(n);
    I{q}.label='True Filter';
    display(I{q}.label)
    [I{q}.n I{q}.P] = FOOPSI2_59(GG,Phat{q},Tim);
    [I{q+numel(qs)}.n I{q+numel(qs)}.P]   = WienerFilt1_2(F,Sim.dt,Phat{q});
    I{q+2*numel(qs)}.n = I{q+numel(qs)}.n; I{q+2*numel(qs)}.n(I{q+2*numel(qs)}.n<0)=0;
    I{q+3*numel(qs)}.n = [diff(F); 0];

    %% compute statistics

%     [I{q}.roc I{q}.auc] = roc2(I{q}.n, D{q}.n);
%     [I{q+numel(qs)}.roc I{q+numel(qs)}.auc] = roc2(I{q+numel(qs)}.n, D{q}.n);
%     xc=corrcoef(I{1}.n,D{1}.n);
%     cc(q,1)=xc(2);
%     xc=corrcoef(I{q+numel(qs)}.n,D{1}.n);
%     cc(q,2)=xc(2);
%     I{q+2*numel(qs)}.n=I{q+numel(qs)}.n; I{q+2*numel(qs)}.n(I{q+2*numel(qs)}.n<0)=0;
%     xc=corrcoef(I{q+2*numel(qs)}.n,D{1}.n);
%     cc(q,3)=xc(2);

    khat0(q,:)=stim'\n;
    khat1(q,:)=stim'\I{1}.n;
    khat2(q,:)=stim'\I{q+numel(qs)}.n;
    khat3(q,:)=stim'\I{q+2*numel(qs)}.n;
    khat4(q,:)=stim'\I{q+3*numel(qs)}.n;

end

%%
mink=min([P.lam(:); khat0(:); khat1(:); khat2(:); khat3(:); khat4(:)]);
maxk=max([P.lam(:); khat0(:); khat1(:); khat2(:); khat3(:); khat4(:)]);

save('kernel.mat')
%% end) plot results
clear Pl
nrows   = 2+Nc;                                  % set number of rows
ncols   = 2;
h       = zeros(nrows,1);
Pl.xlims= [5 Sim.T-101];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 4;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = Sim.Nc;
fnum = figure(1); clf,
Pl.interp = 'latex';

figure(3), clf, hold on
plot(P.lam,'k','LineWidth',Pl.lw);

errorbar(mean(khat0),std(khat0),'--','Color',Pl.gray,'LineWidth',Pl.lw);
errorbar(mean(khat1),std(khat1),'Color',Pl.c(1,:),'LineWidth',Pl.lw)
errorbar(mean(khat2),std(khat2),'-.','Color',Pl.c(2,:),'LineWidth',Pl.lw);
errorbar(mean(khat3),std(khat3),'-.','Color',Pl.c(3,:),'LineWidth',Pl.lw);
errorbar(mean(khat4),std(khat4),'-.','Color',Pl.c(4,:),'LineWidth',Pl.lw);

axis([.9 5.1 mink maxk])
set(gca,'YTick',linspace(mink,maxk,5),'YTickLabel',0:5)
set(gca,'XTick',[1:dim])
xlabel('Kernel Dimension','FontSize',Pl.fs)
ylab=ylabel('Strength','Interpreter','none','FontSize',Pl.fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')


% print fig
wh=[7 5];   %width and height
DirName = '../../figs/';
FileName = 'kernel';
PrintFig(wh,DirName,FileName);
