% this script compares Wiener and Fast Filters
% 
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nNoisy Simulation Fig\n')

% 1) set simulation metadata
Sim.T       = 400;          % # of time steps
Sim.dt      = 0.005;        % time step size
Sim.Plot    = 1;            % whether to plot
Sim.MaxIter = 100;          % max number of iterations
Sim.Mstep   = 1;

% 2) initialize parameters
P.lam   = 8;                % rate, ie, lam*dt=# spikes per second                   
P.tau   = 0.3;              % calcium decay time constant (sec)
P.sig   = 0.3;              % standard deviation of noise (nM)

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T,1));   % simulate spike train  
C = filter(1,[1 -(1-Sim.dt/P.tau)],n);      % calcium concentratin
F = C+P.sig*randn(Sim.T,1);            % fluorescence

% 4) estimate params from real spikes and plot F,C,n
Phat = FastParams1_3(F,C,n,Sim.T,Sim.dt);
P, Phat, figure(10), clf, plot(C), hold on, plot(F,'k'), bar(n)

%% 5) infer spikes and estimate parameters

P2 = P;
P2.lam = P.lam;
% P2.tau = .5*P.tau;
% P2.mu  = 1*P.mu;
P2.sig = P.sig;
Algs=[200 523];                                       % which algorithms within DataComp to use
for m=Algs
    Sim.Alg = m; 
    tic
    I{m}    = DataComp15(F,P2,Sim);
    toc
end
%% 6) plot results

fig=figure(1); clf,
nrows = 3+numel(Algs);
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD
Sim.nticks = 5;
Pl = PlotParams(Sim);

% plot fluorescence data
i=1; subplot(nrows,1,i)
Plot_F(Pl,F)

% plot calcium
i=i+1; subplot(nrows,1,i)
Plot_C(Pl,C)

% plot spike train
i=i+1; subplot(nrows,1,i)
Plot_n(Pl,n)
 
% plot inferred spike trains
for m=Algs
    i=i+1; subplot(nrows,1,i), 
    if m==2
        Pl.label = Pl.WienerLabel;
    elseif m>3 || m>14
        Pl.label = Pl.FastLabel;
    end
    Plot_n_MAP(Pl,I{m}.n/max(I{m}.n))
end

subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)

% print fig
wh=[7 3];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','schem')