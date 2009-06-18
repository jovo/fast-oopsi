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
tau     = 0.5; %round(100*rand(Sim.Nc,1))/100+0.05;   % decay time constant for each cell
P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
dim     = 5;
stim    = rand(Sim.T,dim)'*2;
lam   = [.01 .05 .1 .5 1]; %sin(linspace(0,pi,dim))';


% 3) simulate data
for tt=1:5
    disp(['tt', num2str(tt)])
    clear I
    for q=1:length(lam)
        disp(['lam ', num2str(q)])
        n       = zeros(Sim.T,Sim.Nc);
        C       = n;
        n(1)    = C_0;
        P.lam   = lam(q);
        n(2:end)= poissrnd(P.lam*ones(Sim.T-1,1));    % simulate spike train
        n(n>1)  = 1;
        C       = filter(1,[1 -P.gam],n);               % calcium concentration
        F       = C*P.a' + P.b'+P.sig*randn(Sim.T,Npixs);               % fluorescence

        D{tt,q}.n=n; D{tt,q}.C=C; D{tt,q}.F=F;

        %% infer spikes
        GG=D{q}.F; Tim=Sim;
        Phat{q}=P;
        I{q}.label='True Filter';
        display(I{q}.label)
        I{1,q}.n = FOOPSI_v3_05_01(F',Phat{q},Tim);
        I{2,q}.n = WienerFilt1_2(F,Sim.dt,P);
        I{3,q}.n = I{2,q}.n; I{3}.n(I{3,q}.n<0)=0;
        I{4,q}.n = [diff(F); 0];

        for i=1:4
            snr{q}(tt,i)=mean(I{i,q}.n(D{tt,q}.n==1).^2)/mean(I{i,q}.n(D{tt,q}.n==0).^2);
        end
    end

end

%%
for q=1:length(lam)
    mean_snr(q,:)=mean(snr{q});
    std_snr(q,:) =std(snr{q});
end

save('../../data/snr.mat')

%% end) plot results
load('../../data/snr.mat')
clear Pl
clc
nrows=1; ncols=1;
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
Pl.fs = 18;

subplot(nrows,ncols,nrows*ncols)
% lb=max(mean_snr'*0.99,mean_snr'-std_snr');
% h=errorbar(repmat([1:5],4,1),mean_snr',mean_snr'/100,std_snr');
h=errorbar(mean_snr',std_snr');
set(gca,'YScale','log') %,'YTick',10.^(-5:10),'XTickLabel',[]);
set(h,'LineWidth',2)
ymax=max(mean_snr(:)+std_snr(:));
ymin=max(10^-1,min(mean_snr(:)-lb(:)));
axis([.9 5.1 ymin ymax*1.1])
% set(gca,'YTick',logspace(log(ymin),log(ymax),5),'YTickLabel',0:5)
set(gca,'XTick',1:5,'XTickLabel',lam,'FontSize',Pl.fs)
set(gca,'YTick',10.^(0:5),'YTickLabel',10.^(0:5))
xlabel('$\lambda$','FontSize',Pl.fs,'Interpreter','latex')
ylab=ylabel([{'SNR'}],'Interpreter','none','FontSize',Pl.fs,'Color','k');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')


% print fig
wh=[7 5];   %width and height
DirName = '../../figs/';
FileName = 'snr_sim';
PrintFig(wh,DirName,FileName);
