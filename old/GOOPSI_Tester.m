function I = GOOPSI_Tester(varargin)
% this function takes in either no input (in which case it simulates data),
% or 1 input (which is the fluorescence time series. in either case, it
% runs the GOOPSI filter to infer stuff
fprintf('\nGOOPSI_Tester\n')

% set simulation metadata
Sim.T       = 50;                                   % # of time steps
Sim.dt      = 1/50;                                 % time step size
Sim.freq    = 5;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 20;                                   % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep    = 0;                                   % do M-step
Sim.C_params = 0;                                   % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = 0;                                   % whether to estimate rate governing parameters {b,k}
Sim.h_params = 0;                                   % whether to estimate spike history parameters {h}
Sim.F_params = 0;                                   % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter  = 0;                                   % max # of EM iterartions

% initialize particle filter parameters
P.rate      = 10/(Sim.T*Sim.dt);                    % expected spike rate
P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
P.tau_c     = 2;                                    % calcium decay time constant (sec)
P.A         = 10;                                   % jump size (\mu M)
P.C_0       = 0.1;                                  % baseline [Ca++]
P.C_init    = P.C_0;                                % initial [Ca++]
P.sigma_c   = 1;                                    % standard deviation of noise (\mu M)
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.a         = Sim.dt/P.tau_c;                       % for coding brevity only
P.gamma     = 2e-5;                                 % scaled variance
P.zeta      = 5*P.gamma;                            % constant variance

if Sim.M==1                                         % if there are spike history terms
    P.omega = -0.5;                                 % weight
    P.tau_h = 0.015;                                % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

if nargin==1
    F=z1(varargin{1});
    Sim.T = length(F);
    Sim.dt = 0.031;
else % simulate spikes
    n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
    C=P.C_init*ones(1,Sim.T);
    epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
    for t=2:Sim.T                                           %recursively update calcium
        C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);
    end
    S=Hill_v1(P,C);
    F=(P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(1,Sim.T));
    F(F<=0)=eps;
end

% infer spikes and estimate parameters
Algs = [7];
for m=Algs;
    Sim.Alg = m;
    I{m}    = DataComp14(F,P,Sim);
    if Sim.Alg==7, I{m}.n=I{m}.M.nbar; end
end

%% plot results
fig=figure(1); clf,
row = 1+2*(1-nargin);
nrows = row+length(Algs); 
subplot(nrows,1,1), plot(F), axis('tight'), ylabel('F')
if nargin==0
    subplot(nrows,1,2), plot(C), axis('tight'), ylabel('C')
    subplot(nrows,1,3), bar(n), axis('tight'), ylabel('n')
end
for m=Algs
    row=row+1;
    subplot(nrows,1,row), bar(I{m}.n), axis('tight'), ylabel(I{m}.name)
end