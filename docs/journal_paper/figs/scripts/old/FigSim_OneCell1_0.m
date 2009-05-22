clear all, clc

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size

% 2) initialize spatial filter
mu      = [0 0];
Sigma   = diag([1,1])*1;
Nrows   = 30;
x1      = linspace(-5,5,Nrows); 
x2      = x1;
[X1,X2] = meshgrid(x1,x2);
P.a     = mvnpdf([X1(:) X2(:)],mu,Sigma);
figure(1), clf, imagesc(reshape(P.a,length(x2),length(x1)));

% 3) initialize other params
P.b     = 0*P.a;
P.sig   = .05;                                  % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = 10;                                   % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T-1,1));     % simulate spike train
% n(n>1)=1;
n = [C_0; n];                                   % set initial calcium
C = filter(1,[1 -P.gam],n);                     % calcium concentration
Z = 0*C;
F = C*P.a'+Z*P.b'+P.sig*randn(Sim.T,length(P.a));             % fluorescence

for i=1:Sim.T
    imwrite(reshape(F(i,:),Nrows,Nrows)',['Sim1Cell.tif'],'tif','Compression','none','WriteMode','append')
end

% 4) infer spike train
n=FOOPSI2_56(F,P,Sim);

% 5) plot results