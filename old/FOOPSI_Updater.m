clear, clf, clc, fprintf('\nFOOPSI_Updater\n')

% set simulation metadata
Sim.T       = 100;                                 % # of time steps
Sim.dt      = 0.05;                                % time step size
Sim.Mstep   = true;                                % do M-step
Sim.Plot    = true;                               % plot results as they are obtained
Sim.MaxIter = 10;                                  % max # of EM iterartions

% initialize FOOPSI parameters
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 0.1;                                      % standard deviation of noise (\mu M)
rate    = 10;                                       % expected spike rate
P.lam   = Sim.T/(rate)*Sim.dt;                      % expected jump size ber time bin
P.k     = log(-log(1-rate/Sim.T)/Sim.dt);           % linear filter

% simulate spikes
n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
C=filter(1,[1 -(1-Sim.dt/P.tau)],n);
F=C+P.sig.*randn(Sim.T,1);
F(F<=0)=eps;

% infer spikes and estimate parameters
Algs = [4 8];
for m=Algs;
    Sim.Alg = m;
    P.tau=P.tau/2;
    I{m}    = DataComp14(F,P,Sim);
end

%% plot results
fig=figure(1); clf,
nrows = 3+length(Algs); i=3;

subplot(nrows,1,1), plot(F), axis('tight'), ylabel('F')
subplot(nrows,1,2), plot(C), axis('tight'), ylabel('C')
subplot(nrows,1,3), bar(n), axis('tight'), ylabel('n')
for m=Algs
    i=i+1;
    subplot(nrows,1,i), bar(I{m}.n), axis('tight'), ylabel(I{m}.name)
end