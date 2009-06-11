mu=pi;
T=1e3;

er = exprnd(mu*ones(T,1));
figure(10), clf, hist(er)
muhat = expfit(er)

muhat2 = poissfit(er)