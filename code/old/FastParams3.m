function P = FastParams3(F,C,n,T,dt,gam,con)
% this function estimates the parameters of the model:
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma*C_{t-1} + nu + rho*n_t,   n_t ~ Poisson(n_t; p_t)
%
% if 'Q' is an input into the code, then B is NOT estimated,
% because they are unidentifiable


X  = [gam*C(1:end-1)+n(2:end), 1+0*n(2:end)];              % this is just [gam*C+n, 1]
Y  = F(2:end);
options = optimset('Display','off');
B = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[inf inf],[],options);

if con==0       % unconstrained
    %     B = X\Y;
    P.alpha = B(1);
    P.beta  = B(2);
elseif con==1    % constrained
    P.alpha = B(1);
    P.beta  = B(2);
end

DD      = norm(Y-X*B)^2;                            % difference vector
P.sig   = sqrt(DD/T);

P.lam   = sum(n)/(T*dt);

c       = 1/(2*P.sig^2);                    % scale of variance
P.gam   = gam;
P.l     = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end