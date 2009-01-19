function l = GetLik(F,n,T,dt,P)
    C   = filter(1,[1 -P.gam],n);               % calcium concentratin
    l   = 0.5*T*log(2*pi*P.sig^2) + norm(F-P.a*C-P.b)^2/(2*P.sig^2)...
        - T*log(P.lam*dt) + P.lam*dt*sum(n);
end