% this script does makes a figure showing how the smc-em approach
% outperforms can incorporate prior information to refine the spike train
% estimates by doing the following:
%
% 1) load some data
% 2) set simulation metadata (eg, dt, T, # particles, etc.)
% 3) initialize parameters
% 4) infers spikes using a variety of approaches
% 5) plots results
clear, clc, close all, fprintf('\nStimulus Fig\n')

%% 1) get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load '/Users/joshyv/Research/data/rafa/brendon/last_2_weeks/080218a/D080218a.mat';
j=24;                                               % experiment number

D{j}.spt    = Get_spt(D{j}.V);                      % get spike times
D{j}.n      = SubSampleSpikes(D{j},1);              % make spike train

%% 2) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T       = min(length(D{j}.n),length(D{j}.F));   % # of time steps
Sim.dt      = D{j}.dt_o;                            % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 1;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.n       = D{j}.n; Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

Sim.Mstep   = true;                                 % do M-step
Sim.C_params= true;                                % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params= true;                                % whether to estimate spike history parameters {h}
Sim.F_params= false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 50;                                   % max # of EM iterartions

%% 3) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize barrier and wiener filter parameters
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 1.65;                                        % calcium decay time constant (sec)
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.tau_c     = 1.6472;
P.A         = 20.8036;
P.C_0       = 19.2702;                                   % baseline [Ca++]
P.C_init    = 0.02;                                 % initial [Ca++]
P.sigma_c   = 8.8704;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 2;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 2e-5;                                 % scaled variance
P.zeta      = 5*P.gamma;                                 % constant variance

if Sim.M==1                                         % if there are spike history terms
    P.omega = -0.0941;                                   % weight
    P.tau_h = 0.015;                                    % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

%% 4) infer spikes and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Algs=[4 5];                                           % which algorithms within DataComp to use
Sim.freq    = 1;
tvec_o      = Sim.freq:Sim.freq:Sim.T;              % observation time steps
for m=Algs
    Sim.Alg = m;
    if any(m==[2 4]),
        F       = D{j}.F(tvec_o);
        Tim     = Sim;
        Tim.dt  = Sim.dt*Sim.freq;
        Tim.T   = Sim.T/Sim.freq;
        P.lam   = Sim.T/(sum(D{j}.n(1:Sim.T))*P.A)*Sim.dt;  % expected jump size ber time bin
    elseif m==7
        F           = D{j}.F(1:Sim.T);
        Sim.StimDim = 1;                                    % # of stimulus dimensions
        Sim.x       = ones(Sim.StimDim,Sim.T);
        Tim         = Sim;
        P.k         = 1.1065; %[log(-log(1-sum(D{j}.n)/Sim.T)/Sim.dt)];% linear filter
    elseif m==9
        F           = D{j}.F(1:Sim.T);
        Sim.StimDim = 2;                                    % # of stimulus dimensions
        x           = SubSampleStim(D{j},1)';               % incorporate stimulus
        Sim.x       = ones(Sim.StimDim,numel(x));
        Sim.maxx    = max(x);
        Sim.minx    = min(x);
        Sim.x(2,:)  = x/Sim.maxx;                     % normalize stimuls
        Tim         = Sim;
        P.k         = [log(-log(1-sum(D{j}.n)/Sim.T)/Sim.dt); 1];% linear filter
    elseif any(m==[5 6])
        P.A=0.1; P.sig=10;
        F       = D{j}.F(tvec_o);
        Tim     = Sim;
        Tim.dt  = Sim.dt*Sim.freq;
        Tim.T   = Sim.T/Sim.freq;
        P.lam   = Sim.T/(sum(D{j}.n(1:Sim.T))*P.A)*Sim.dt;  % expected jump size ber time bin
    end
    I{m} = DataComp13(F,P,Tim);
end

%% 5) plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig=figure(2); clf,
nrows=numel(Algs)+1;
gray  = [.75 .75 .75];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                     % interpreter for axis labels
xlims = [80 Sim.T-3];               % xmin and xmax for current plot
fs=12;                              % font size
ms=15;                              % marker size for real spike
sw=1.7;                             % spike width
sw2=1.7;                              % spike width for wiener
lw=2;                               % line width
I{2}.name=[{'Optimal'}; {'Linear Filter'}];
I{4}.name=[{'Optimal'}; {'Nonnegative Filter'}];
I{5}.name=[{'fPPR'}; {'Filter'}];
I{7}.name=[{'Optimal Nonlinear'}; {'Particle Filter'}];
I{9}.name=[{'Generalized Optimal'}; {'Nonlinear Filter'}];
Sim.n=D{j}.n;
spt=find(Sim.n==1);
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(D{j}.FrameStartTime>i,1);
end

% plot real fluorescence
i=1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(D{j}.F(tvec_o))+1,'-k','LineWidth',2,'MarkerSize',ms*.8);
stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Real Fluorescence'}; {'and Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 2.1])

% plot real spike train
% i=i+1; subplot(nrows,1,i), hold on
% stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
% ylab=ylabel([{'Real'}; {'Spikes'}],'Interpreter',inter,'FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% set(gca,'YTick',0:2,'YTickLabel',[])
% set(gca,'XTick',secs,'XTickLabel',[])
% axis([xlims 0 1])

% % plot stimulus
% if any(Algs==9)
%     i=i+1; subplot(nrows,1,i), hold on
%     plot(Sim.x(2,:),'Color','k','LineWidth',lw)
%     ylab=ylabel([{'Stimulus'};{['(\mu' 'A)']}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
%     set(gca,'YTick',[-.5 .5],'YTickLabel',[{'-36'}; {'+35'}]); %round([Sim.minx/2 Sim.maxx/2]*1e-3))
%     set(gca,'XTick',secs,'XTickLabel',[])
%     axis([xlims min(Sim.x(2,:)) max(Sim.x(2,:))])
% end

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if any(m==[7 9])
        stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        axis([xlims 0 1])
    elseif m==2
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(zeros(size(Sim.n)),'k')
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate
        
                
        neg = find(n_est<=0);
        stem(tvec_o(neg),n_est(neg),'Marker','none','LineWidth',sw2,'Color',col(1,:))        
        pos = find(n_est>0);
        stem(tvec_o(pos),n_est(pos),'Marker','none','LineWidth',sw2,'Color',col(2,:))

%         stem(Sim.freq:Sim.freq:Sim.T,n_est,'Marker','none','LineWidth',8,'Color',col(1,:));
        %         area(tvec_o,n_est,'FaceColor',col(1,:),'EdgeColor',col(1,:))
        axis([xlims min(n_est) max(n_est)])
    elseif any(m==[4 5 6]) 
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        plot(zeros(size(Sim.n)),'k')
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate
        stem(tvec_o,n_est,'Marker','none','LineWidth',sw2,'Color',col(2,:))
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
%     plot(C,'Color',col(2,:),'LineWidth',2)
%     set(gca,'YTick',1,'YTickLabel',[])
%     ylab=ylabel([{'Inferred'};{'Calcium'}],'Interpreter',inter,'FontSize',fs);
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
%     axis([xlims min(C) max(C)+.1])
% end

set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','NIPSdata')
