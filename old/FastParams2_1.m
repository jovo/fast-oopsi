function P = FastParams2_1(F,C,n,T,dt)


P.lam  = sum(n)/(T*dt);

W   = [C(1:end-1) 1+0*C(1:end-1)];              % this is just [C 1]
Y   = F(2:end)-n(2:end);
a = W\Y;
if a(1)>1, a(1)=1; elseif a(1)<0, a(1)=0; end           % make sure 'a' is within bounds

P.a     = a(1);
P.tau   = dt/(1-a(1));
P.mu    = a(2);

DD   = (F-C-P.mu)'*(F-C-P.mu);                            % difference vector
P.sig= sqrt(DD/T);

c  = 1/(2*P.sig^2);                    % scale of variance
P.lik= 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end