function P = FastParams2_7(F,C,n,T,dt,Q)
% this function estimates the parameters of the model:
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma*C_{t-1} + nu + rho*n_t,   n_t ~ Poisson(n_t; p_t)
%
% if 'Q' is an input into the code, then B is NOT estimated,
% because they are unidentifiable

if nargin==5
    X  = [C(1:end-1) 1+0*n(2:end) n(2:end)];              % this is just [C 1 n]
    Y  = F(2:end);                                      
    [B,fval,exitflag,output,lambda] = quadprog(X'*X,-X'*Y,[],[],[],[],[0 -inf -inf],[1 inf inf]);
    P.gam   = B(1);
    P.nu    = B(2);
    P.rho   = B(3);
elseif nargin==6
    X       = [C 1+0*n n];              % this is just [C 1 n]
    Y       = F;                                      
    [B,fval,exitflag,output,lambda] = quadprog(X'*X,-X'*Y,[],[],[],[],[0 -inf -inf],[1 inf inf]);
    P.gam   = B(1);
    P.nu    = B(2);
    P.rho   = B(3);
    P.nu    = mean(Y-Q.gam*C-Q.rho*n);
%     X       = [C 1+0*n n];              % this is just [C 1 n]
%     Y       = F(2:end);                                      
%     P.gam   = Q.gam;
%     P.rho   = Q.rho;
%     P.nu    = Q.nu;
%     B       = [P.gam P.nu P.rho]';
end

DD          = norm(Y-X*B)^2;                            % difference vector
P.sig       = sqrt(DD/T);

P.lam  = sum(n)/(T*dt);

c           = 1/(2*P.sig^2);                    % scale of variance
P.lik       = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end