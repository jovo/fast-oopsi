function l = GetNonlinLik(F,C,n,T,dt,P)
%     C   = filter(1,[1 -P.gam],n);               % calcium concentration
    l   = 0.5*T*log(2*pi*P.sig^2) + norm(F-P.a*C./(C+P.k_d)-P.b)^2/(2*P.sig^2)...
        - T*log(P.lam*dt) + P.lam*dt*sum(n);
end