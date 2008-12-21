function P = FastParams2_3(F,C,n,T,dt,Q)
% this function estimates the parameters of the model:
% n_t ~ Poisson(n_t; p_t)
% C_t = a C_{t-1} + n_t, where a=(1-dt/tau), and A is fixed at 1.
% F_t = C_t + eps_t, eps_t ~ N(mu,sig^2)
% 
% if 'Q' is an input into the code, then both a and mu are NOT estimated,
% but 'Q' indicates that {n_t} is estimated as well, and therefore these
% params become unidentifiable

P.lam  = sum(n)/(T*dt);

if nargin==5
    W       = [C(1:end-1) 1+0*C(1:end-1)];              % this is just [C 1]
    Y       = F(2:end)-n(2:end);
    a       = W\Y;
    if a(1)>1, a(1)=1; elseif a(1)<0, a(1)=0; end           % make sure 'a' is within bounds
    P.a     = a(1);
    P.tau   = dt/(1-a(1));
    P.mu    = a(2);
elseif nargin==6
    P.tau   = Q.tau;
    P.mu    = Q.mu; %sum(F-C)/T;
end

DD          = (F-C-P.mu)'*(F-C-P.mu);                            % difference vector
P.sig       = sqrt(DD/T);

c           = 1/(2*P.sig^2);                    % scale of variance
P.lik       = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end