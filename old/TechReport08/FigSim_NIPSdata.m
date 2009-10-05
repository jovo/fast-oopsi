% this script does makes a figure showing how the smc-em approach
% outperforms linear algorithms with an example data set by doing:
%
% 1) load some data and get some parameters using real data
% 2) set simulation metadata (eg, dt, T, # particles, etc.)
% 3) initialize parameters
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nSaturation Fig\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load '/Users/joshyv/Research/data/rafa/brendon/last_2_weeks/080218a/D080218a.mat';
j=14;                                               % experiment number

% arrange stuff for plotting purposes
freq        = 1;
D{j}.spt    = Get_spt(D{j}.V);                      % get spike times
D{j}.n      = SubSampleSpikes(D{j},freq);              % make spike train
Sim.n       = D{j}.n; %Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD
T           = min(length(D{j}.n),length(D{j}.F));   % # of time steps
F           = NaN*ones(T,1);
tvec_o      = freq:freq:T;
F(tvec_o)   = D{j}.F(1:T);
figure(1), clf, plot(F,'.'), hold on, stem(Sim.n,'k')

% [x C F] = GetDataParams(D{j}); %x0 = [A a C_0 k_d F_max];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T_o     = T;                                    % # of time steps
Sim.dt      = D{j}.dt_o;                            % time step size
Sim.freq    = freq;                                    % # of time steps between observations
Sim.T       = Sim.T_o*Sim.freq;                     % # of observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 100;                                  % # of particles
Sim.M       = 1;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep   = true;                                 % do M-step
Sim.C_params = true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                                % whether to estimate spike history parameters {h}
Sim.F_params = false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 30;                                    % max # of EM iterartions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize barrier and wiener filter parameters
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 2;                                        % calcium decay time constant (sec)
P.lam   = Sim.T/(sum(D{j}.n(1:Sim.T))*P.A)*Sim.dt;  % expected jump size ber time bin
P.sig   = 0.1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-sum(D{j}.n)/Sim.T)/Sim.dt);% linear filter
P.tau_c     = 1.6;%Sim.dt/x(2);
P.A         = 50;
P.C_0       = 18.5;                                   % baseline [Ca++]
P.C_init    = 0.02;                                 % initial [Ca++]
P.sigma_c   = 9;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 2;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 1e-5;                                 % scaled variance
P.zeta      = 5e-5;                                 % constant variance

if Sim.M==1                                         % if there are spike history terms
    P.omega = -0.5;                                   % weight
    P.tau_h = 0.015;                                    % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) infer spikes and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Algs=[2 4];                                       % which algorithms within DataComp to use
Sim.freq = 1;
for m=Algs
    Sim.Alg = m;
    if m==9,
        P.k=[P.k; 0];
        Sim.StimDim = 2;                                    % # of stimulus dimensions
        x           = SubSampleStim(D{j},1)';               % incorporate stimulus
        Sim.x       = ones(Sim.StimDim,numel(x));
        Sim.maxx    = max(x);
        Sim.minx    = min(x);
        Sim.x(2,:)  = x/Sim.maxx;                     % normalize stimuls
    end
    I{m}    = DataComp13(D{j}.F(1:Sim.T),P,Sim);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig=figure(2); clf,
if any(Algs==7) || any(Algs==9), nrows = numel(Algs)+2;
else nrows = numel(Algs)+2; end
gray  = [.75 .75 .75];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'none';                     % interpreter for axis labels
xlims = [1 Sim.T];              % xmin and xmax for current plot
fs=12;                              % font size
ms=15;                              % marker size for real spike
sw=1.8;                             % spike width
lw=2;                               % line width
I{2}.name=[{'Optimal'}; {'Linear Filter'}];
I{7}.name=[{'Optimal'}; {'Nonlinear Filter'}];
spt=find(D{j}.n>=1); Sim.n=zeros(size(Sim.n)); Sim.n(spt)=spt;
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-2);
for i=1:Nsec-2
    secs(i) = find(D{j}.FrameStartTime>i,1);
end

% plot real data
i=1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(D{j}.F(tvec_o)),'-k','LineWidth',2,'MarkerSize',ms);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Real'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1.1*max(D{j}.F(1:xlims(2)))])

% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(D{j}.n,'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Spike'}; {'Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
axis([xlims 0 max(D{j}.n)])

% plot stimulus
if any(Algs==9)
    i=i+1;
    subplot(nrows,1,i), hold on
    plot(Sim.x(2,:),'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'};{'(mA)'}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[-.5 .5],'YTickLabel',round([Sim.minx/2 Sim.maxx/2]*1e-3))
    axis([xlims min(Sim.x(2,:)) max(Sim.x(2,:))])
end

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if any(m==[7 9])
        stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        axis([xlims 0 1])
    else
        stem(spt,D{j}.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        n_est   = I{m}.n; n_est   = max(D{j}.n)*n_est/max(n_est);   %normalize estimate

        neg = find(n_est<=0);
        stem(neg,n_est(neg),'Marker','none','LineWidth',sw,'Color',col(1,:))
        pos = find(n_est>0);
        stem(pos,n_est(pos),'Marker','none','LineWidth',sw,'Color',col(2,:))

        axis([xlims min(n_est) max(n_est)])
    end

    hold off,
    ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',secs,'XTickLabel',[])
end

% % plot calcium
% if any(Algs==7) && ~any(Algs==9)
%     subplot(nrows,1,i+1), hold on
%     C = I{m}.M.Cbar/I{m}.P.k_d;
%     hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(1,:) I{m}.M.Cptiles(end,Sim.T:-1:1)]/I{m}.P.k_d,ccol(2,:));
%     set(hfill,'edgecolor',ccol(2,:))
%     plot(C(1:end-1),'Color',col(2,:),'LineWidth',2)
%     set(gca,'YTick',1,'YTickLabel',[])
%     ylab=ylabel([{'Inferred'};{'Calcium'}],'Interpreter',inter,'FontSize',fs);
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
%     axis([xlims min(C) 1.1*max(C(1:xlims(2)))])
% end


set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 4];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','NIPSdata')