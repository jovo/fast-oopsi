function P = FastParams2_4(F,C,n,T,dt,Q)
% this function estimates the parameters of the model:
% n_t ~ Poisson(n_t; p_t)
% C_t = a C_{t-1} + n_t, where a=(1-dt/tau), and A is fixed at 1.
% F_t = C_t + eps_t, eps_t ~ N(mu,sig^2)
% 
% if 'Q' is an input into the code, then both a and mu are NOT estimated,
% but 'Q' indicates that {n_t} is estimated as well, and therefore these
% params become unidentifiable


if nargin==5
    X       = [C(1:end-1) 1+0*n(2:end) n(2:end)];              % this is just [C 1 n]
    Y       = F(2:end);
    A       = X\Y;
    if A(1)>1, A(1)=1; elseif A(1)<0, A(1)=0; end           % make sure 'a' is within bounds
    P.gamma = A(1);
    P.nu    = A(2);
    P.rho   = A(3);
elseif nargin==6
    P.gamma = Q.gamma;
    P.nu    = Q.nu;
    P.rho   = Q.rho;
end

DD          = norm(Y-X*A)^2;                            % difference vector
P.sig       = sqrt(DD/T);

P.lam  = sum(n)/(T*dt);

c           = 1/(2*P.sig^2);                    % scale of variance
P.lik       = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end