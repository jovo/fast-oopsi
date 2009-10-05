function P = FastParams2_5(F,C,n,T,dt,Q)
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
    B       = X\Y;
    if B(1)>1, B(1)=1; elseif B(1)<0, B(1)=0; end           % make sure 'a' is within bounds
    P.tau   = dt/(1-B(1));
    P.C_b   = B(2)/(1-B(1));
    P.A     = B(3)*P.tau/dt;
    P.gamma = B(1);
    P.nu    = B(2);
    P.rho   = B(3);
    
    H=X'*X;
    f=-X'*Y;
    B = quadprog(H,f,[],[],[],[],[0 -inf -inf],[1 inf inf]);

    P.tau2  = dt/(1-B(1));
    P.C_b2  = B(2)/(1-B(1));
    P.A2    = B(3)*P.tau2/dt;

    P.gamma2= B(1);
    P.nu2   = B(2);
    P.rho2  = B(3);

elseif nargin==6
    P.tau   = Q.tau;
    P.C_b   = Q.C_b;
    P.A     = Q.A;
    P.gamma = Q.gamma;
    P.nu    = Q.nu;
    P.rho   = Q.rho;
end

DD          = norm(Y-X*B)^2;                            % difference vector
P.sig       = sqrt(DD/T);

P.lam  = sum(n)/(T*dt);

c           = 1/(2*P.sig^2);                    % scale of variance
P.lik       = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end