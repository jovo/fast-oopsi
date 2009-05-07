% this script compares Wiener and Fast Filters
%
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nNonlinear model\n')

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 0;                               % max number of iterations

% 2) initialize parameters
P.a     = .1;
P.b     = 0;
P.sig   = .001;                                  % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = ((linspace(0,100,Sim.T-1)).*(sin(linspace(0,12*pi,Sim.T-1))+1))';                                   % rate-ish, ie, lam*dt=# spikes per second
P.k_d   = 180;
P.rho   = 50;

% 3) simulate data
% n = poissrnd(P.lam*Sim.dt);     % simulate spike train
% n = poissrnd(P.lam*Sim.dt*ones(Sim.T-1,1));     % simulate spike train
spt = [100 200:10:220 300:5:320 400:2:420];
n = zeros(Sim.T-1,1);
n(spt) = 1;
% n(n>1)=1;
% n(n>1)=1;
n = [C_0; n];                                   % set initial calcium
C = filter(1,[1 -P.gam],n*P.rho);                     % calcium concentration
S = C./(C+P.k_d);
F = P.a*S+P.b+P.sig*randn(Sim.T,1);             % fluorescence


X=0:.01:max(C); P.n=1;
figure(10), clf, nrows=4; i=0;
i=i+1; subplot(nrows,1,i), plot(F,'k'), axis('tight')
hold on, plot(P.a*S+P.b), axis('tight')
i=i+1; subplot(nrows,1,i), hold on, plot(z1(P.a*S+P.b)), axis('tight')
plot(z1(C),'r'), axis('tight')
i=i+1; subplot(nrows,1,i), bar(n), axis('tight')
i=i+1; subplot(nrows,1,i), semilogx(Hill_v1(P,X));

% 4) estimate params from real spikes
P.case=2;
Phat = FastParams3_3(F,C,n,Sim.T,Sim.dt,P);
display(Phat);

%% 5) infer spikes and estimate parameters
P2=P;
P2.lam=Phat.lam;

% initialize parameters
for q=1:2
    if q==1
        [I{q}.n I{q}.P] = FOOPSI2_55(F,P2,Sim);
        I{q}.P.label = 'Linear';
    elseif q==2
        Sim.n=I{1}.n;
        [I{q}.n I{q}.P] = NFOOPSI1_0(F,P2,I{1}.P,Sim);
        I{q}.P.label = 'Nonlinear';
    end
end

%% 6) plot results
fig     = figure(1); clf,
nrows   = 4+q;                                  % set number of rows
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

% plot inferred calcium
i=i+1; h(2) = subplot(nrows,1,i);
Pl.label = [{'Inferred'}; {'Calcium'}];
Pl.color = Pl.gray;
plot(C,'k'), hold on
plot(max(C)*z1(filter(1,[1 -P.gam],I{1}.n*P.rho)),'b')
plot(filter(1,[1 -P.gam],I{2}.n*P.rho),'r')
axis('tight')
% X=0:.01:max(C); P.n=1;
% i=i+1; subplot(nrows,1,i), semilogx(Hill_v1(P,X));


subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
% linkaxes(h,'x')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','nonlin')