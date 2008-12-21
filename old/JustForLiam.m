% function I = OOPSI_Tester(varargin)
% this function takes in either no input (in which case it simulates data),
% or 1 input (which is the fluorescence time series. in either case, it
% runs some OOPSI filters to infer stuff
fprintf('\nOOPSI_Tester\n')

% set simulation metadata
Sim.T       = 100;
Sim.dt      = 0.031;                                % median frame duration
Sim.MaxIter = 30;                                  % max # of EM iterartions
Sim.Plot    = 1;                                   % plot results with each iteration

% initialize parameters
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 0.1;                                      % standard deviation of noise (\mu M)
Nsp     = 50;                                       % expected spike rate
P.lam   = Sim.T/(Nsp*Sim.dt);                       % expected jump size ber time bin

% simulate
n=rand(Sim.T,1)<P.lam*Sim.dt*exp(-P.lam*Sim.dt);
C=filter(1,[1 -(1-Sim.dt/P.tau)],n);
F=C+P.sig.*randn(Sim.T,1);

% infer spikes and estimate parameters
[I.n I.P]   = FOOPSI_v1_4(F,P,Sim);

%% plot results
figure(1); clf,
subplot(411), plot(F), axis('tight'), ylabel('F')
subplot(412), plot(C), axis('tight'), ylabel('C')
subplot(413), bar(n), axis('tight'), ylabel('n')
subplot(414), bar(I.n), axis('tight'), ylabel('filter')
