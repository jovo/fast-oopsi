function I = OOPSI_Tester(varargin)
% this function takes in either no input (in which case it simulates data),
% or 1 input (which is the fluorescence time series. in either case, it
% runs some OOPSI filters to infer stuff
clc, fprintf('\nOOPSI_Tester\n')

% set simulation meta-data
if nargin==1                                        % if input fluorescence vector
    F       = z1(varargin{1});                      % normalize F to be between 0 and 1
    Sim.T   = length(F);                            % # of time steps
else
    Sim.T   = 1000;
end
Sim.dt      = 0.050;                                % median frame duration

% M-step meta-data 
Sim.Mstep    = 1;                                   % do M-step
Sim.MaxIter  = 20;                                % max # of EM iterartions
Sim.Plot     = 1;                                   % plot results with each iteration

% initialize fast parameters
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 0.1;                                      % standard deviation of noise (\mu M)
Nsp     = Sim.T/5;                                       % expected spike rate
P.lam   = Sim.T/(Sim.dt*Nsp);                       % expected jump size ber time bin

% particle filtering meta-data
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.N       = 99;                                   % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(Sim.StimDim,Sim.T);              % stimulus

% particle filtering M-step meta-data
Sim.C_params = 0;                                   % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = 0;                                   % whether to estimate rate governing parameters {b,k}
Sim.h_params = 0;                                   % whether to estimate spike history parameters {h}
Sim.F_params = 0;                                   % whether to estimate observation parameters {alpha,beta,gamma,zeta}

% initialize particle filter parameters
P.rate      = Nsp/(Sim.T*Sim.dt);                   % expected spike rate
P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
P.tau_c     = 2;                                    % calcium decay time constant (sec)
P.A         = 10;                                   % jump size (\mu M)
P.C_0       = 0.1;                                  % baseline [Ca++]
P.C_init    = P.C_0;                                % initial [Ca++]
P.sigma_c   = 1;                                    % standard deviation of noise (\mu M)
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 2;                                    % F_max
P.beta      = 0;                                    % F_min
P.a         = Sim.dt/P.tau_c;                       % for coding brevity only
P.gamma     = 1e-5;                                 % scaled variance
P.zeta      = 5*P.gamma;                            % constant variance

if Sim.M==1                                         % if there are spike history terms
    P.omega = -0.5;                                 % weight
    P.tau_h = 0.015;                                % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

if nargin==0                                                % if no input fluorescence vector    
    n=rand(Sim.T,1)<P.lam*Sim.dt*exp(-P.lam*Sim.dt);
    C=filter(1,[1 -(1-Sim.dt/P.tau)],n);
    F=C+P.sig.*randn(Sim.T,1);
%     F(F<=0)=eps;
%     n=rand(Sim.T,1)<1-exp(-exp(P.k*Sim.x')*Sim.dt);         % sample spikes
%     C=P.C_init*ones(Sim.T,1);                               % extize calcium vector
%     epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      % generate noise on calcium
%     for t=2:Sim.T                                           % recursively update calcium
%         C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);
%     end
%     S=Hill_v1(P,C);                                         % apply nonlinearity for fluorescence
%     F=(P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(Sim.T,1)); % add noise
%     F(F<=0)=eps;                                            % make sure F_t >=0 for all t
end

% infer spikes and estimate parameters
Algs = [4 41];
for m=Algs;
    Sim.Alg = m; 
%     if m==41, 
%         P.tau=P.tau*2; 
%         P.sig=P.sig*2; 
%         P.lam=200;
%         Sim.Mstep=1; 
%     end
    I{m}    = DataComp14(F,P,Sim);
end

%% plot results
figure(1); clf,
row = 1+2*(1-nargin);
nAlgs = numel(Algs);
nrows = row+nAlgs; 
subplot(nrows,1,1), plot(F), axis('tight'), ylabel('F')
if nargin==0
    subplot(nrows,1,2), plot(C), axis('tight'), ylabel('C')
    subplot(nrows,1,3), bar(n), axis('tight'), ylabel('n')
end
for m=Algs
    row=row+1;
    subplot(nrows,1,row), bar(I{m}.n), axis('tight'), ylabel(I{m}.name)
end