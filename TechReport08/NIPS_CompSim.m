% this script does makes a figure showing how the smc-em approach
% outperforms other approaches even for just noisy data, by doing:
%
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, close all, fprintf('\nExample Simulation Fig\n')

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T       = 405;                                  % # of time steps
Sim.dt      = 0.025;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep   = false;                                % do M-step
Sim.C_params = true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                                % whether to estimate spike history parameters {h}
Sim.F_params = false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 0;                                    % max # of EM iterartions

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize barrier and wiener filter parameters

rate    = .025;
spt     = poissrnd(rate*ones(Sim.T,1)); %poissrnd(linspace(0,10,Sim.T));
P.rate  = sum(spt)/(Sim.T*Sim.dt);                % expected spike rate
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 0.5;                                        % calcium decay time constant (sec)
P.lam   = 1/(P.rate*P.A);                           % expected jump size ber time bin
P.sig   = 0.4;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
P.tau_c     = P.tau;
P.A         = 1;
P.C_0       = 0;                                   % baseline [Ca++]
P.C_init    = 0;                                % initial [Ca++]
P.sigma_c   = P.sig;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 0e-5;                                 % scaled variance
P.zeta      = 5e-5;                            % constant variance
P.a         = Sim.dt/P.tau_c;

%% 3) simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
n=spt;
C=P.C_init*ones(Sim.T,1);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
for t=2:Sim.T                                           %recursively update calcium
    C(t)  = (1-P.a)*C(t-1) + P.a*P.C_0 + P.A*n(t);% + epsilon_c(t);
end
R.C=C;
eps_t=P.sig*randn(Sim.T,1);
F=C+eps_t;
% S=Hill_v1(P,R.C);
% F=P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*eps_t;
% F(F<=0)=eps;
% figure(2), clf,
% subplot(311), plot(z1(F)+1); hold on, stem(n),
% subplot(312), plot(S+sqrt(P.gamma*S+P.zeta).*eps_t);
% subplot(313), plot(C*max(S)/6+sqrt(P.zeta).*eps_t);
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

%% 4) infer spikes and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Algs=[2 4 5];                                       % which algorithms within DataComp to use
Sim.freq = 1;
for m=Algs
    Sim.Alg = m;
    I{m}    = DataComp13(F,P,Sim);
end

%% 5) plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig=figure(1); clf,
ncols = 2;
nrows = numel(Algs)+2;
gray  = [.75 .75 .75];              % define gray color
col   = [1 0 0; 0 .5 0];            % define colors for mean
ccol  = col+.8; ccol(ccol>1)=1;     % define colors for std
inter = 'none';                     % interpreter for axis labels
xlims = [2 Sim.T-2];              % xmin and xmax for current plot
fs=16;                              % font size
ms=20;                              % marker size for real spike
sw=2.5;                             % spike width
lw=2;                               % line width
% make xticks
Nsec = floor(Sim.T*Sim.dt);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(Sim.tvec>=i,1);
end
I{2}.name=[{'Optimal'}; {'Linear Filter'}];
I{4}.name=[{'Optimal'}; {'Nonnegative Filter'}];
I{5}.name=[{'fPPR'}; {'Filter'}];
I{7}.name=[{'Inferred'}; {'Spikes'}];

minn=min(n);
maxn=max(n);
for m=Algs
    minn=min(minn,min(I{m}.n));
    maxn=max(maxn,max(I{m}.n));
end
% plot fluorescence data
i=1; subplot(nrows,ncols,i), hold on
plot(z1(F),'k','LineWidth',2);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[0 1],'YTickLabel',[0 1])
set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
axis([xlims 0 1.0])
title(['Rate = ' num2str(rate/Sim.dt) ' Hz'],'FontSize',fs+4)

% plot real spike train
i=i+2; subplot(nrows,ncols,i), hold on
% plot(spt,'k','linewidth',lw)
stem(spt,'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'# Spikes'}; {'per Frame'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[0:max(spt)],'FontSize',fs)
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 max(spt)])

% plot inferred spike trains
for m=Algs
    i=i+2;
    subplot(nrows,ncols,i), hold on,

    if any(m==[7 9])
        %         subplot(nrows,ncols,i), hold on,
        stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        axis([xlims 0 1])
    else
        %         subplot(nrows,ncols,i), hold on,

        %         plot(n,'k','LineWidth',lw)

        %         neg = find(I{m}.n<=0);
        %         nneg = 0*ones(size(I{m}.n));
        %         nneg(neg)=I{m}.n(neg);
        %         plot(nneg,'-','Color',col(1,:),'LineWidth',lw)
        %
        %         pos = find(I{m}.n>0);
        %         npos = 0*ones(size(I{m}.n));
        %         npos(pos)=I{m}.n(pos);
        %         plot(npos,'-','Color',col(2,:),'LineWidth',lw)
        %
        %         plot(zeros(size(I{m}.n)),'k','linewidth',lw)
        %         axis('tight')
        
        sptt=spt; sptt(spt==0)=[];
        pos=find(spt>0);
        stem(pos,spt(pos),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')

        if m==2
        neg = find(I{m}.n<=0);
        stem(neg,I{m}.n(neg),'Marker','none','LineWidth',sw,'Color',col(1,:))
        end
        pos = find(I{m}.n>0);
        stem(pos,I{m}.n(pos),'Marker','none','LineWidth',sw,'Color',col(2,:))

        plot(zeros(size(I{m}.n)),'k','linewidth',lw)
        axis([xlims minn maxn])
    end

    hold off,
    ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[0 1 2],'YTickLabel',[],'FontSize',fs) %round([0 (maxn-minn)/2+minn maxn]*10)/10
    set(gca,'XTick',secs,'XTickLabel',[])
    set(gca,'XTickLabel',[])
    box off
end

subplot(nrows,ncols,nrows*ncols-1)
set(gca,'XTick',secs,'XTickLabel',1:Nsec,'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

%% 6) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T       = 105;                                  % # of time steps
Sim.dt      = 0.025;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep   = false;                                % do M-step
Sim.C_params = true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                                % whether to estimate spike history parameters {h}
Sim.F_params = false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 0;                                    % max # of EM iterartions

%% 7) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize barrier and wiener filter parameters

rate    = 5;
spt     = poissrnd(rate*ones(Sim.T,1)); %poissrnd(linspace(0,10,Sim.T));
% spt     = conv2(spt,exp(-(Sim.T/2-[0:Sim.T]).^2/10)','same');
% spt     = round(rate*spt/max(spt));
P.rate  = numel(spt)/(Sim.T*Sim.dt);                % expected spike rate
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 0.5;                                        % calcium decay time constant (sec)
P.lam   = 1/sqrt(P.rate*P.A);                           % expected jump size ber time bin
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
P.tau_c     = P.tau;
P.A         = 1;
P.C_0       = 1;                                   % baseline [Ca++]
P.C_init    = 100;                                % initial [Ca++]
P.sigma_c   = P.sig;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 0e-5;                                 % scaled variance
P.zeta      = 5e-5;                            % constant variance
P.a         = Sim.dt/P.tau_c;

%% 8) simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
n=spt;
C=P.C_init*ones(Sim.T,1);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
for t=2:Sim.T                                           %recursively update calcium
    C(t)  = (1-P.a)*C(t-1) + P.a*P.C_0 + P.A*n(t);% + epsilon_c(t);
end
R.C=C;
eps_t=P.sig*randn(Sim.T,1);
F=C+eps_t;
% S=Hill_v1(P,R.C);
% F=P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*eps_t;
% F(F<=0)=eps;
% figure(2), clf,
% subplot(311), plot(z1(F)+1); hold on, stem(n),
% subplot(312), plot(S+sqrt(P.gamma*S+P.zeta).*eps_t);
% subplot(313), plot(C*max(S)/6+sqrt(P.zeta).*eps_t);
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

%% 9) infer spikes and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.freq = 1;
for m=Algs
    Sim.Alg = m;
    I{m}    = DataComp13(F,P,Sim);
end

%% 10) plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig=figure(1); clf,
% ncols = 2;
% nrows = 4;
% gray  = [.75 .75 .75];              % define gray color
% col   = [1 0 0; 0 .5 0];            % define colors for mean
% ccol  = col+.8; ccol(ccol>1)=1;     % define colors for std
% inter = 'none';                     % interpreter for axis labels
% fs=16;                              % font size
% ms=20;                              % marker size for real spike
% sw=2.5;                             % spike width
% lw=2;                               % line width
figure(1)
xlims = [2 Sim.T-2];              % xmin and xmax for current plot
% make xticks
% Nsec = floor(Sim.T*Sim.dt*10);
% secs = zeros(1,Nsec-1);
% for i=1:Nsec
%     secs(i) = find(Sim.tvec>=i,1);
% end
secs=[20:20:Sim.T];
% I{2}.name=[{'Optimal'}; {'Linear Filter'}];
% I{4}.name=[{'Optimal'}; {'Nonnegative Filter'}];
% I{5}.name=[{'fPPR'}; {'Filter'}];
% I{7}.name=[{'Inferred'}; {'Spikes'}];

minn=min(n);
maxn=max(n);
for m=Algs
    minn=min(minn,min(I{m}.n));
    maxn=max(maxn,max(I{m}.n));
end
% plot fluorescence data
i=2; subplot(nrows,ncols,i), hold on
plot(z1(F),'k','LineWidth',2);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
% ylab=ylabel([{'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[0 1],'YTickLabel',[0 1])
set(gca,'XTick',secs,'XTickLabel',[],'FontSize',fs)
axis([xlims 0 1.0])
title(['Rate = ' num2str(rate/Sim.dt) ' Hz'],'FontSize',fs+4)

% plot real spike train
i=i+2; subplot(nrows,ncols,i), hold on
plot(spt,'k','linewidth',lw)
% stem(spt,'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
% ylab=ylabel([{'Real'}; {'Spikes'}],'Interpreter',inter,'FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[0 rate rate*2],'FontSize',fs)
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 rate*2])

% plot inferred spike trains
for m=Algs
    i=i+2;
    subplot(nrows,ncols,i), hold on,

    if any(m==[7 9])
        %         subplot(nrows,ncols,i), hold on,
        stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        axis([xlims 0 1])
    else
        %         subplot(nrows,ncols,i), hold on,
        plot(n,'k','LineWidth',lw)

        neg = find(I{m}.n<=0);
        nneg = 0*ones(size(I{m}.n));
        nneg(neg)=I{m}.n(neg);
        plot(nneg,'-','Color',col(1,:),'LineWidth',lw)

        pos = find(I{m}.n>0);
        npos = 0*ones(size(I{m}.n));
        npos(pos)=I{m}.n(pos);
        plot(npos,'-','Color',col(2,:),'LineWidth',lw)

        plot(zeros(size(I{m}.n)),'k','linewidth',lw)

        %         stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        %
        %         neg = find(n_est<=0);
        %         stem(neg,n_est(neg),'Marker','none','LineWidth',sw,'Color',col(1,:))
        %
        %         pos = find(n_est>0);
        %         stem(pos,n_est(pos),'Marker','none','LineWidth',sw,'Color',col(2,:))
        %         axis('tight')
        axis([xlims 0 rate*2])
    end

    hold off,
    %     ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    %     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[0 rate rate*2],'YTickLabel',[],'FontSize',fs)
    set(gca,'XTick',secs,'XTickLabel',[])
    set(gca,'XTickLabel',[])
    box off
end

subplot(nrows,ncols,nrows*ncols)
set(gca,'XTick',secs,'XTickLabel',[{'0.2'}; {'0.4'}; {'0.6'}; {'0.8'}; {'1'}],'FontSize',fs)
xlabel('Time (msec)','FontSize',fs)

% print fig
wh=[11 7];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print('-depsc','NIPScomp')