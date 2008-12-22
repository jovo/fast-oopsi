% this script compares Wiener and Fast Filters
% 
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results
% 
% differs from FigSim_Schem1 in that we estimate beta in 
% F_t = alpha C_t + beta + sigma*epsilon_t
% whereas in FigSim_Schem1 we assumed we knew it.

clear, clc, fprintf('\nNoisy Simulation Fig\n')

% 1) set simulation metadata
Sim.T       = 500;          % # of time steps
Sim.dt      = 0.005;        % time step size
Sim.Plot    = 1;            % whether to plot
Sim.MaxIter = 100;          % max number of iterations

% 2) initialize parameters
P.alpha = 1;
P.beta  = 0;
P.sig   = 0.01;
tau     = 0.8;
C_b     = 0;
P.gamma = 1-Sim.dt/tau;         % tau       = dt/(1-gamma)
P.nu    = Sim.dt*C_b/tau;       % baseline  = nu/(1-gamma)
P.rho   = 0.2;                  % Sim.dt*A/tau; A = rho*tau/dt
P.lam   = 4;                    % rate, ie, lam*dt=# spikes per second                   

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T,1));       % simulate spike train  
C = filter(P.rho,[1 -P.gamma],n+P.nu/P.rho);    % calcium concentratin
F = P.alpha*C+P.beta+P.sig*randn(Sim.T,1);      % fluorescence

Phat = FastParams2_6(F,C,n,Sim.T,Sim.dt);
P, Phat, figure(10), clf, plot(P.alpha*C+P.beta), hold on, plot(F,'k'), bar(n)

%% 5) infer spikes and estimate parameters

P2 = P;
P2.lam = 2*P.lam;
% P2.tau = .5*P.tau;
% P2.mu  = 0*P.mu;
% P2.sig = P.sig;
Algs=5.25; %[2.11 5.24];                                       % which algorithms within DataComp to use
q   = 1;
for m=Algs
    Sim.Alg = m; 
    tic
    I{q}    = DataComp16(F,P2,Sim);
    toc
    q=q+1;
end
%% 6) plot results

fig=figure(1); clf,
nrows = 3+numel(Algs);
Pl.xlims=[101 500];
Pl.nticks=5;
Pl.n=double(n); Pl.n(Pl.n==0)=NaN;
% Sim.n = double(n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD
% Sim.nticks = 5;
Pl = PlotParams(Pl);

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
    q=1;
for m=Algs
    i=i+1; subplot(nrows,1,i), 
    if m>=2 && m <300
        Pl.label = Pl.WienerLabel;
    elseif m>=500 || m<600
        Pl.label = Pl.FastLabel;
    end
    Plot_n_MAP(Pl,I{q}.n/max(I{q}.n))
    q=q+1;
end

subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',(Pl.XTicks-Pl.XTicks(1))*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','schem')