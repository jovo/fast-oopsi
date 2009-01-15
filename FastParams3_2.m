function P = FastParams3_2(F,C,n,T,dt,Q)
% this function estimates the parameters of the model:
% F_t = a*C_t + b + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma*C_{t-1} + nu + rho*n_t,   n_t ~ Poisson(n_t; p_t)
%
% Inputs 
%   F:  fluorescence time series
%   C:  calcium
%   n:  spike train
%   T:  total number of time steps
%   dt: time step size
%   Q:  initial parametersm, and which ones should we estimate
% 
% Outputs
%   P:  updated parameters

options = optimset('Display','off');
P       = Q;                                                % initialize all parameters
switch Q.case
    case 0 % do nothing
        P.name = 'true';
    case 1 % estimate a and b, unconstrained
        X   = [C 1+0*C];
        Y   = F;
        B   = X\Y;
        P.a = B(1);
        P.b = B(2);
        P.name = 'unconstrained';
    case 2 % estimate a, b > 0
        X   = [C 1+0*C];
        Y   = F;
        B   = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[inf inf],[Q.a, Q.b],options);
        P.a = B(1);
        P.b = B(2);     
        P.name = 'a,b>0';
    case 3 % estimate a and b, 0<a<1, b>0
        X   = [C 1+0*C];
        Y   = F;
        B   = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[1 inf],[],options);
        P.a = B(1);
        P.b = B(2);      
        P.name = '0<a<1, b>0';
    case 4 % estimate only a, unconstrained
        X   = C;
        Y   = F;
        B   = X\Y;
        P.a = B(1);
        P.name = 'a unconstrained';
    case 5 % a>0
        X   = C;
        Y   = F;
        B   = quadprog(X'*X,-X'*Y,[],[],[],[],0,inf,[],options);
        P.a = B(1);
        P.name = 'a>0';
    case 6 % 0<a<1
        X   = C;
        Y   = F;
        B   = quadprog(X'*X,-X'*Y,[],[],[],[],0,1,[],options);
        P.a = B(1);
        P.name = '0<a<1';
    case 7 % estimate b > 0
        X   = F-P.a*C;
        Y   = 0*F;
        P.b = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[inf inf],[],options);
        P.name = 'b>0';
end

DD      = norm(F-P.a*C-P.b)^2;
P.sig   = sqrt(DD/T);
P.lam   = sum(n)/(T*dt);

c       = 1/(2*P.sig^2);                    % scale of variance
P.l     = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end