function P = FastParams2_8(F,C,n,T,dt,Q)
% this function estimates the parameters of the model:
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma*C_{t-1} + nu + rho*n_t,   n_t ~ Poisson(n_t; p_t)
%
% if 'Q' is an input into the code, then B is NOT estimated,
% because they are unidentifiable

X  = [C(1:end-1) 1+0*n(2:end) n(2:end)];              % this is just [C 1 n]
Y  = F(2:end);
if nargin==5
    B = quadprog(X'*X,-X'*Y,[],[],[],[],[0 -inf -inf],[1 inf inf]);
    P.gam   = B(1);
    P.nu    = B(2);
    P.rho   = B(3);
elseif nargin==6
% XX  = [C 1+0*n n];              % this is just [C 1 n]
% YY  = F;
% 
%     options = optimset('Display','off');
%     B = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0 0],[1 inf inf],[],[],options);
%     P.gam   = B(1);
%     P.nu    = B(2);
%     P.rho   = B(3);
    %     YY = [F(1:T-1)-Q.rho*n(1:T-1) - Q.nu];
    %     XX = C(2:T);
    %     P.gam   = quadprog(XX'*XX,-XX'*YY,[],[],[],[],[0 -inf -inf],[1 inf inf]);
    %     P.nu    = Q.nu;
    %     P.rho   = Q.rho;
    %     B       = [P.gam P.nu P.rho]';
    P.gam   = Q.gam;
    P.rho   = Q.rho;
    P.nu    = Q.nu;
    B       = [P.gam P.nu P.rho]';
end

DD      = norm(Y-X*B)^2;                            % difference vector
P.sig   = sqrt(DD/T);

P.lam   = sum(n)/(T*dt);

c       = 1/(2*P.sig^2);                    % scale of variance
P.l     = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end