P.a     = 1;                % scale
P.b     = 0;                % bias
tau     = 0.5;              % decay time constant for each cell
P.gam   = 1-V.dt./tau;      % set gam
lams    = [1 3 5 10 30 50 100];                % rate
P.sig   = 0.2;
mse.f   = zeros(1,length(lams));
mse.w   = mse.f;
V.T     = 300;
V.dt    = 1/30;
for j=1:length(lams)
    P.lam=lams(j);
    for i=1%:10 % trials
        V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));   % simulate spike train
        V.C     = filter(1,[1 -P.gam(1)],V.n);         % calcium concentration
        V.F     = P.a*V.C+P.b+P.sig*randn(V.T,1);
        fast    = fast_oopsi(V.F,V,P);
        fast    = fast(2:end-1);
        fast    = fast/max(fast)*max(V.n);
        wiener  = wiener_oopsi(V.F,V.dt,P);
        wiener  = wiener(2:end-1);
        wiener  = wiener/max(abs(wiener))*max(V.n);
        mse.f(i,j) = sum((fast-V.n(2:end-1)).^2);
        mse.w(i,j) = sum((wiener-V.n(2:end-1)).^2);
        wt      = wiener;
        wt(wt<0)= 0;
        mse.wt(i,j) = sum((wt-V.n(2:end-1)).^2);
    end
    figure(2), clf,
    plot(z1(V.F)+max(V.n)),
    hold all,
    bar(V.n);
    stem(fast)
    plot(wiener)
    [P.lam mse.f(1,j) mse.w(1,j) mse.wt(1,j)]
    pause 
end
