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
% 2_3: don't know, something else buggy i think
% 2_4: reparameterized a bit.

clear, clc, fprintf('\ndemo\n')

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 30;                               % max number of iterations

% 2) initialize parameters
P.alpha = 1;
P.beta  = 0;
P.sig   = .3;                                   % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;                         
P.lam   = 10;                                   % rate-ish, ie, lam*dt=# spikes per second                   

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T-1,1));     % simulate spike train  
n = [C_0; n];                                   % set initial calcium
C = filter(1,[1 -P.gam],n);                     % calcium concentration
F = P.alpha*C+P.beta+P.sig*randn(Sim.T,1);      % fluorescence

% 4) estimate params from real spikes and plot F,C,n
Phat = FastParams3_1(F,C,n,Sim.T,Sim.dt);
orderfields(P), orderfields(Phat), figure(10), clf, 
plot(Phat.alpha*filter(1,[1 -Phat.gam],n)+ Phat.beta,'linewidth',2), 
hold on, plot(F,'k'), bar(-n), axis([1000 1500 -1 5])
Sim.dt/(1-Phat.gam), tau
%% 5) infer spikes and estimate parameters

% initialize parameters
P2 = P;
P2.lam  = 2*P.lam;
P2.sig  = 2*P.sig;
% P2.alpha= 2*P.alpha;
% P2.beta = 2*P.beta;

I.name      = [{'Fast'}; {'Filter'}];
[I.n I.P]   = FOOPSI2_51(F,P,Sim);

%% 6) plot results

orderfields(P), orderfields(Phat), orderfields(I.P)

fig=figure(1); clf,
nrows = 4;                                      % set number of rows
Pl.xlims=[1 Sim.T];                             % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n=double(n); Pl.n(Pl.n==0)=NaN;              % store spike train for plotting
Pl = PlotParams(Pl);                            % generate a number of other parameters for plotting

% plot fluorescence data
i=1; h(1) = subplot(nrows,1,i);
Pl.label = 'Fluorescence';          
Pl.color = 'k';
Plot_X(Pl,F);

% plot calcium
i=i+1; h(2) = subplot(nrows,1,i);
Pl.label = 'Calcium';
Pl.color = Pl.gray;
Plot_X(Pl,C);

% plot spike train
i=i+1; h(3) = subplot(nrows,1,i);
maxn=max(n(Pl.xlims(1):Pl.xlims(2)));
Plot_n(Pl,n/maxn);
 
% plot inferred spike trains
i=i+1; h(4) = subplot(nrows,1,i);
Pl.label = Pl.FastLabel;
maxn=max(I.n(Pl.xlims(1):Pl.xlims(2)));
Plot_n_MAP(Pl,I.n);

% xlabel properly
subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',(Pl.XTicks-Pl.xlims(1))*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','schem')