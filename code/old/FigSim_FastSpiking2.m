% this script compares Wiener and Fast Filters
% 
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results
% 
% 2_0: we estimate beta in F_t = alpha C_t + beta + sigma*epsilon_t
% 2_1: totally reparameterized. see section 3
% 2_2: allow for initial condition on C

clear, clc, fprintf('\nNoisy Simulation Fig\n')

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.05;                            % time step size
Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 0;                               % max number of iterations

% 2) initialize parameters
P.sig   = 0.3;                                  % stan dev of noise
tau     = 0.8;                                  % decay time constant
C_b     = 0;                                    % baseline
P.gam   = 1-Sim.dt/tau;                         % tau       = dt/(1-gamma)
P.nu    = Sim.dt*C_b/tau;                       % baseline  = nu/(1-gamma)
P.rho   = 1;                                  % jump size
P.C_0   = 0;                                    % calcium at first time step
P.lam   = 100;                                    % rate-ish, ie, lam*dt=# spikes per second                   

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*(sin(linspace(-4*pi, 8*pi,Sim.T-1))'+2));     % simulate spike train  
n = [P.C_0/P.rho; n];                           % set initial calcium
C = filter(P.rho,[1 -P.gam],n+P.nu/P.rho);      % calcium concentration
F = C+P.sig*randn(Sim.T,1);                     % fluorescence

% 4) estimate params from real spikes and plot F,C,n
figure(10), clf, plot(C), hold on, plot(F,'k'), bar(n)
Phat = FastParams2_7(F,C,n,Sim.T,Sim.dt);
P, Phat, 

%% 5) infer spikes and estimate parameters

% initialize parameters
P2 = P;
P2.lam = 2*P.lam;
P2.sig = 2*P.sig;

Algs=[2.12 5.243];                              % set which algorithms within DataComp to use
q=1;                                            % count algs
for m=Algs
    Sim.Alg = m;                                % set alg
    tic
    I{q}    = DataComp16(F,P2,Sim);             % actually use algorithm 'm' to infer spike train and learn parameters
    toc
    q=1+1;                                      % update counter
end

%% 6) plot results
fig=figure(1); clf,
nrows = 3+numel(Algs);                          % set number of rows
Pl.xlims=[100 Sim.T-2];                           % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n=double(n); Pl.n(Pl.n==0)=NaN;              % store spike train for plotting
Pl = PlotParams(Pl);                            % generate a number of other parameters for plotting

% plot fluorescence data
i=1; subplot(nrows,1,i)
Pl.label = 'Fluorescence';          
Pl.color = 'k';
Plot_X(Pl,F)

% plot calcium
i=i+1; subplot(nrows,1,i)
Pl.label = 'Calcium';
Pl.color = Pl.gray;
Plot_X(Pl,C)

% plot spike train
i=i+1; subplot(nrows,1,i)
maxn=max([max(n) max(I{1}.n) max(I{2}.n)]);
Plot_n(Pl,n/maxn)
 
% plot inferred spike trains
q=0;
for m=Algs
    q=q+1; i=i+1; subplot(nrows,1,i),            
    if m>=2 && m <3
        Pl.label = Pl.WienerLabel;
    elseif m>=5 || m<6
        Pl.label = Pl.FastLabel;
    end
    Plot_n_MAP(Pl,I{q}.n/maxn)
end

subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',(Pl.XTicks-Pl.xlims(1))*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','FastSpiking')