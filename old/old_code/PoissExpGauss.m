clear
T=1e5;
dt=.005;
lam=100;

n=poissrnd(lam*dt*ones(T,1));

nhat(1) = poissfit(n);
nhat(2) = expfit(n);
[mu si] = normfit(n);
% nhat(3) = mu;
hist(n,100)
[lam*dt nhat mu si^2]
