% function I = OOPSI_Tester(varargin)
% this function takes in either no input (in which case it simulates data),
% or 1 input (which is the fluorescence time series. in either case, it
% runs some OOPSI filters to infer stuff
clc, fprintf('\nOOPSI_Tester\n')

% set simulation metadata
if nargin==1                                        % if input fluorescence vector
    F       = z1(varargin{1});                      % normalize F to be between 0 and 1
    Sim.T   = length(F);                            % # of time steps
else
    Sim.T   = 500;
end
Sim.dt      = 0.031;                                % median frame duration
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.N       = 99;                                   % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(Sim.StimDim,Sim.T);              % stimulus

Sim.Mstep    = 1;                                   % do M-step
Sim.C_params = 0;                                   % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = 0;                                   % whether to estimate rate governing parameters {b,k}
Sim.h_params = 0;                                   % whether to estimate spike history parameters {h}
Sim.F_params = 0;                                   % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter  = 100;                                  % max # of EM iterartions
Sim.Plot     = 1;                                   % plot results with each iteration

% initialize FOOPSI parameters
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 0.1;                                      % standard deviation of noise (\mu M)
Nsp     = 100;                                       % expected spike rate
% P.lam   = -1/Sim.dt*lambertw(-Nsp/Sim.T);
P.lam   = Nsp/(Sim.T*Sim.dt);                       % expected jump size ber time bin

% % initialize particle filter parameters
% P.rate      = Nsp/(Sim.T*Sim.dt);                   % expected spike rate
% P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
% P.tau_c     = 2;                                    % calcium decay time constant (sec)
% P.A         = 10;                                   % jump size (\mu M)
% P.C_0       = 0.1;                                  % baseline [Ca++]
% P.C_init    = P.C_0;                                % initial [Ca++]
% P.sigma_c   = 1;                                    % standard deviation of noise (\mu M)
% P.n         = 1.0;                                  % hill equation exponent
% P.k_d       = 200;                                  % hill coefficient
% P.alpha     = 2;                                    % F_max
% P.beta      = 0;                                    % F_min
% P.a         = Sim.dt/P.tau_c;                       % for coding brevity only
% P.gamma     = 1e-5;                                 % scaled variance
% P.zeta      = 5*P.gamma;                            % constant variance
% 
% if Sim.M==1                                         % if there are spike history terms
%     P.omega = -0.5;                                 % weight
%     P.tau_h = 0.015;                                % time constant
%     P.sigma_h = 0.01;                               % stan dev of noise
% end

% if nargin==0                                        % if no input fluorescence vector
    n=rand(Sim.T,1)<P.lam*Sim.dt;%*exp(-P.lam*Sim.dt);
    C=filter(1,[1 -(1-Sim.dt/P.tau)],n);
    F=C+P.sig.*randn(Sim.T,1);
    %     n=rand(Sim.T,1)<1-exp(-exp(P.k*Sim.x')*Sim.dt);         % sample spikes
    %     C=P.C_init*ones(Sim.T,1);                               % extize calcium vector
    %     epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      % generate noise on calcium
    %     for t=2:Sim.T                                           % recursively update calcium
    %         C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);
    %     end
    %     S=Hill_v1(P,C);                                         % apply nonlinearity for fluorescence
    %     F=(P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(Sim.T,1)); % add noise
    %     F(F<=0)=eps;                                            % make sure F_t >=0 for all t
% end

% infer spikes and estimate parameters
Algs = [17];
for m=Algs;
    Sim.Alg = m;
%     E.lam   = P.lam*2;
    E.sig   = P.sig*2;
    E.tau   = P.tau/1.5;
    I{m}    = DataComp14(F,E,Sim);
end

%% plot results
figure(1); clf,
row = 1+2;%*(1-nargin);
nAlgs = numel(Algs);
nrows = row+nAlgs;
ax(1)=subplot(nrows,1,1); plot(z1(F)), axis('tight'),
ylab=ylabel('F');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% if nargin==0
    ax(2)=subplot(nrows,1,2); plot(z1(C)), axis('tight'), ylab=ylabel('Calcium'); 
    set(gca,'YTickLabel',[])
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    ax(3)=subplot(nrows,1,3); bar(n), axis('tight'), ylab=ylabel('n');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% end
for m=Algs
    row=row+1;
    ax(row)=subplot(nrows,1,row); bar(I{m}.n), axis('tight'), ylab=ylabel(I{m}.name);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
end
linkaxes(ax,'xy')