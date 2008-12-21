function P = FastParams1(F,C,n,T,dt)

P.lam  = sum(n)/(T*dt);

W       = C(1:end-1);              
Y       = F(2:end)-n(2:end);
P.a    = W\Y;
if P.a(1)>1, P.a(1)=1; elseif P.a(1)<0, P.a(1)=0; end           % make sure 'a' is within bounds
P.tau= dt/(1-P.a(1));

DD   = (F-C)'*(F-C);                            % difference vector
P.sig= sqrt(DD/T);

P.c  = 1/(2*P.sig^2);                    % scale of variance
P.lik= 0.5*T*log(2*pi*P.sig^2) + P.c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood
P.mu = 0;

end