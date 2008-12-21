function P = FastParams3_1(F,C,n,T,dt,Q)
% this function estimates the parameters of the model:
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma*C_{t-1} + nu + rho*n_t,   n_t ~ Poisson(n_t; p_t)
%
% if 'Q' is an input into the code, then B is NOT estimated,
% because they are unidentifiable

if nargin==5
    X  = [C(1:end-1) 1+0*n(2:end)];              % this is just [C, 1]
    Y  = F(2:end)-n(2:end);
    options = optimset('Display','off');
    B = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[1 inf],[],options);
    P.gam   = B(1);
    P.beta  = B(2);
    DD      = norm(Y-X*B)^2;                            % difference vector
    P.alpha = 1;
    P.lam   = sum(n)/(T*dt);
else
    X  = Q.gam*C(1:end-1)+n(2:end);              % this is just [C, 1]
    Y  = F(2:end);
    options = optimset('Display','off');
    P.alpha = quadprog(X'*X,-X'*Y,[],[],[],[],[0],[inf],[],options);
    %     P.gam   = Q.gam;
    %     P.beta  = mean(F-C);
    %     if P.beta<0, P.beta=0; end
    P.gam   = Q.gam;
    P.beta  = 0;
    DD      = norm(F-P.alpha*C-P.beta)^2;
    %     P.lam   = sum(n(2:T))/(T*dt);
    m       = n(2:end)>.75*max(n);
    P.lam   = sum(m)/(T*dt);
end

P.sig   = sqrt(DD/T);
P.lam   = sum(n)/(T*dt);


c       = 1/(2*P.sig^2);                    % scale of variance
P.l     = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end