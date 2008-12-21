clear, clc, fprintf('\nJust for Adam\n')

load AdamData; %F=F/max(F);
% set simulation metadata
Sim.T       = length(F);                                 % # of time steps
Sim.dt      = 0.031;                                % time step size
Sim.Mstep   = true;                                % do M-step
Sim.Plot    = true;                                % plot results as they are obtained
Sim.MaxIter = 10;                                  % max # of EM iterartions

% initialize FOOPSI parameters
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 30;                                      % standard deviation of noise (\mu M)
rate    = 20;                                       % expected spike rate
P.lam   = Sim.T/(rate)*Sim.dt;  % expected jump size ber time bin
P.k     = log(-log(1-rate/Sim.T)/Sim.dt);           % linear filter

% % simulate spikes
% n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
% C=zeros(Sim.T,1);
% for t=2:Sim.T                                           %recursively update calcium
%     C(t)  = (1-Sim.dt/P.tau)*C(t-1) + n(t);   
% end
% F=C+P.sig.*randn(Sim.T,1);
% F(F<=0)=eps;

% infer spikes and estimate parameters
[n_est P]=FOOPSI_v1_1(F,P,Sim);

%% plot results
fig=figure(1); clf,nrows = 4;
subplot(nrows,1,1), plot(F), axis('tight'), ylabel('F')
% subplot(nrows,1,2), plot(C), axis('tight'), ylabel('C')
% subplot(nrows,1,3), bar(n), axis('tight'), ylabel('n')
subplot(nrows,1,4), bar(n_est), axis('tight'), ylabel('FOOPSI')