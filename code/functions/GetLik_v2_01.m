function post = GetLik_v2_01(DD,n,P,Meta)

sound(3*sin(linspace(0,90*pi,2000)))
lik     = -Meta.T*Meta.Np*log(2*pi*P.sig^2)/2 -1/(2*P.sig^2)*DD;
prior   = Meta.T*sum(P.lam*Meta.dt) - Meta.dt*P.lam'*sum(n)';
post    = lik + prior;