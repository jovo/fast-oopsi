% this script compares Wiener and Fast Filters
%
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results
%
% 2_0: we estimate b in F_t = a C_t + b + sig*epsilon_t
% 2_1: totally reparameterized. see section 3
% 2_2: allow for initial condition on C
% 2_3: don't know, something else buggy i think
% 2_4: reparameterized a bit.

clear, clc, fprintf('\nTrying crap out\n')

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 0;                                % max number of iterations

% 2) initialize parameters
P.a     = .1;
P.b     = .1;
P.sig   = .02;                                    % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = 10;                                   % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T-1,1));     % simulate spike train
n(n>1)=1;
n = [C_0; n];                                   % set initial calcium
C = filter(1,[1 -P.gam],n);                     % calcium concentration
F = P.a*C+P.b+P.sig*randn(Sim.T,1);             % fluorescence

% 4) estimate params from real spikes and plot F,C,n
figure(10), clf,
nn=n/P.a; nn(nn==0)=NaN;
subplot(8,1,1), hold on, plot(F,'k'), bar(-nn), axis([1 Sim.T min(min(F),-P.a) max(F)])

Q = P
for q=1%:7
    Q.case=q-1;
    Phat{q} = FastParams3_2(F,C,n,Sim.T,Sim.dt,Q);
    display(Phat{q});
    subplot(8,1,q+1),
    hold on
    plot(F,'k')
    plot(Phat{q}.a*filter(1,[1 -Phat{q}.gam],n)+ Phat{q}.b,'linewidth',1),
    axis([1 Sim.T -P.a max(F)])
end

%% 5) infer spikes and estimate parameters

% initialize parameters
P2 = P;
% P2.a    = 10;
% P2.b    = 15;
% P2.lam  = 2*P.lam;
% P2.sig  = 2*P.sig;
b=[0:.25*P.b:3*P.b];
reps=length(b)+1;
for q=1:reps
    if q==reps, 
        P2.case=7; 
        P2.b=0;
        Sim.MaxIter=30;
    else
        P2.case=0;
        P2.b=b(q);
    end
    tic
    I{q}.name       = [{'Fast'}; {'Filter'}];
    [I{q}.n I{q}.P] = Bscrew(F,P2,Sim);
    toc
end

%% 6) plot results

nrows   = 3+reps;                                    % set number of rows
h       = zeros(nrows);
Pl.xlims= [5 Sim.T];                           % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;              % store spike train for plotting
Pl      = PlotParams(Pl);                            % generate a number of other parameters for plotting

% plot likelihoods
figure(2), clf
for r=1:reps
    subplot(nrows,1,3+r);
    Pl.label = I{r}.name;
    try
        plot(I{r}.P.l(2:end)), axis('tight')
    catch ME
        display('ass')
    end
end

fig     = figure(1); clf,
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
Plot_n(Pl,n);
% title(['a=',num2str(Phat{1}.a),', b=',num2str(Phat{1}.b), ', sig=',num2str(Phat{1}.sig), ' lam=',num2str(Phat{1}.lam)])

% plot inferred spike trains
for r=1:reps-1
    i=i+1; h(3+r) = subplot(nrows,1,i);
    Pl.label = num2str(I{r}.P.b);
    Plot_n_MAP(Pl,I{r}.n);
%     title();
%     title(['a=',num2str(I{r}.P.a),', b=',num2str(I{r}.P.b), 'sig=',num2str(I{r}.P.sig),...
%         ', lam=',num2str(I{r}.P.lam), ', n=',num2str(max(I{r}.n))])
set(gca,'TickLength',[0.001 0.001])
end


subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
% linkaxes(h,'x')
% linkaxes([h(end-1), h(end)])

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','schem')