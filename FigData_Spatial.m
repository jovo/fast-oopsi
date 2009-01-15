%% this file compares various filters for real data.  specifically D080218a
%% trial 24

load FilteredData_D080218a24

% 1) set simulation metadata
Sim.T       = D.T;                              % # of time steps
Sim.dt      = D.dt_o;                            % time step size

% 2) initialize parameters
P.a     = .1;
P.b     = .1;
P.sig   = .03;                                  % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = 30;                                   % rate-ish, ie, lam*dt=# spikes per second
P.l     = 1e99;                                 % initialize likelihood

% 3) estimate params from real spikes
P.case=2;

for j=1:100
    C = filter(1,[1 -P.gam],D.n);                     % calcium concentration
    P = FastParams3_2(¡F,C,D.n,Sim.T,Sim.dt,P);
end
%%
display(Phat);


Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 30;                               % max number of iterations




% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T-1,1));     % simulate spike train
% n(n>1)=1;
n = [C_0; n];                                   % set initial calcium
C = filter(1,[1 -P.gam],n);                     % calcium concentration
F = P.a*C+P.b+P.sig*randn(Sim.T,1);             % fluorescence

% 4) estimate params from real spikes
P.case=2;
Phat = FastParams3_2(F,C,n,Sim.T,Sim.dt,P);
display(Phat);

%% 5) infer spikes and estimate parameters

% initialize parameters
P2 = P;
P2.a    = P.a/2;
P2.b    = P.b/2;
P2.lam  = 2*P.lam;
P2.sig  = 2*P.sig;

for q=1:2
    P2.case=2;
    I{q}.name       = [{'Fast'}; {'Filter'}];
    if q==1, Sim.thresh=0;
    elseif q==2, Sim.thresh=1; P2.name = [{'Thr'}];
    end
    tic
    [I{q}.n I{q}.P] = FOOPSI2_53(F,P2,Sim);
    toc
end

%% 6) plot results
fig     = figure(1); clf,
nrows   = 3+q;                                    % set number of rows
h       = zeros(nrows);
Pl.xlims= [5 Sim.T];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting

for r=1:q;
    display(max(diff(I{q}.P.l)))
end

% plot likelihoods
figure(2), clf
for r=1:q
    subplot(nrows-3,1,r);
    Pl.label = I{r}.name;
    plot(I{r}.P.l(2:end)), axis('tight')
end

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
title(['a=',num2str(Phat.a),', b=',num2str(Phat.b), ', sig=',num2str(Phat.sig), ' lam=',num2str(Phat.lam)])

% plot inferred spike trains
for r=1:q
    i=i+1; h(3+r) = subplot(nrows,1,i);
    Pl.label = I{r}.P.name;
    Plot_n_MAP(Pl,I{r}.n);
    title(['a=',num2str(I{r}.P.a),', b=',num2str(I{r}.P.b), 'sig=',num2str(I{r}.P.sig),...
        ', lam=',num2str(I{r}.P.lam), ', n=',num2str(max(I{r}.n))])
end


subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
% linkaxes(h,'x')
% linkaxes([h(end-1), h(end)])

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SimThresh')