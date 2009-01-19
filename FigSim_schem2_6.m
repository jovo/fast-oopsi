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
% 2_5: debugged a bit.

clear, clc, fprintf('\nWhy thresholding is useful\n')

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.Plot    = 1;                                % whether to plot
% Sim.MaxIter = 99;                               % max number of iterations

% 2) initialize parameters
P.a     = .1;
P.b     = 0;
P.sig   = .07;                                  % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = 10;                                   % rate-ish, ie, lam*dt=# spikes per second

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
for q=1:9
    if q==1,
        Sim.MaxIter=0;
        I{q}.P.label   = 'no estimate'; 
        P2.label = 'no estimate';    
    else
        Sim.MaxIter=29;
        P2.case = q-2;
        % P2.a    = P.a/2;
        % P2.b    = P.b/2;
        P2.lam  = 2*P.lam;
        P2.sig  = .2*P.sig;
    end
    [I{q}.n I{q}.P] = FOOPSI2_55(F,P2,Sim);
end


%% 6) plot results
for r=1:q
    display(max(diff(I{q}.P.l(2:end))))
end

% plot likelihoods
figure(2), clf
for r=1:q
    subplot(3+q,1,3+r);
    Pl.label = I{r}.P.label;
    plot(I{r}.P.l(2:end),'.-'), axis('tight')
end

% paper fig
fig     = figure(1); clf,
nrows   = 3+q;                                  % set number of rows
h       = zeros(nrows,1);
Pl.xlims= [5 Sim.T];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 2;

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
    Pl.label = I{r}.P.label;
    Plot_n_MAP(Pl,I{r}.n);
    title(['a=',num2str(I{r}.P.a),', b=',num2str(I{r}.P.b), ', sig=',num2str(I{r}.P.sig),...
        ', lam=',num2str(I{r}.P.lam), ', n_{max}=',num2str(max(I{r}.n))])
end


subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
linkaxes(h,'x')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','schem')