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

clear, clc, fprintf('\nNoisy Simulation Fig\n')

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 100;                              % max number of iterations

% 2) initialize parameters

P.alpha = 1;                                    % scale
P.beta  = 0;                                    % offset
P.sig   = 0.1;                                  % stan dev of noise
tau     = 0.8;                                  % decay time constant
C_b     = 0;                                    % baseline
P.gamma = 1-Sim.dt/tau;                         % tau       = dt/(1-gamma)
P.nu    = Sim.dt*C_b/tau;                       % baseline  = nu/(1-gamma)
P.rho   = 0.5;                                  % Sim.dt*A/tau; A = rho*tau/dt
P.lam   = 5;                                    % rate, ie, lam*dt=# spikes per second                   

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T,1));       % simulate spike train  
C = filter(P.rho,[1 -P.gamma],n+P.nu/P.rho);    % calcium concentration
F = P.alpha*C+P.beta+P.sig*randn(Sim.T,1);      % fluorescence

% 4) estimate params from real spikes and plot F,C,n
Phat = FastParams2_6(F,C,n,Sim.T,Sim.dt);
P, Phat, figure(10), clf, plot(C), hold on, plot(F,'k'), bar(n)

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
Pl.xlims=[101 Sim.T];                           % time steps to plot
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
Plot_n(Pl,n)
 
% plot inferred spike trains
q=0;
for m=Algs
    q=q+1; i=i+1; subplot(nrows,1,i),            
    if m>=2 && m <3
        Pl.label = Pl.WienerLabel;
    elseif m>=5 || m<6
        Pl.label = Pl.FastLabel;
    end
    Plot_n_MAP(Pl,I{q}.n/max(I{q}.n))
end

subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',(Pl.XTicks-Pl.XTicks(1))*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','schem')