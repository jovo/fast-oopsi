clear, clc

Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.Nc      = 2;                               % # cells


P.sig   = 0.05;                                  % stan dev of noise
C_0     = 0;
tau     = round(100*rand(Sim.Nc,1))/100+0.05;                                 % decay time constant
P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
P.lam   = round(10*rand(Sim.Nc,1))+5; %10*ones(Sim.Nc,1);                                   % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n=zeros(Sim.T,Sim.Nc);
C=n;
for i=1:Sim.Nc
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(i)*Sim.dt*ones(Sim.T-1,1));    % simulate spike train
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));               % calcium concentration
end
figure(3), ncols=3;
subplot(1,3,1), imagesc(C)
CC=C';
CC=CC(:);
subplot(1,3,2),plot(CC)
Ctemp=reshape(CC,Sim.Nc,Sim.T);
subplot(1,3,3), imagesc(Ctemp')
norm(Ctemp'-C)

