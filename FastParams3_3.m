function [P l] = FastParams3_3(F,C,n,T,dt,P)
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
%   P:  initial parameters, and which ones should we estimate
%
% Outputs
%   P:  updated parameters
%   l:  lik

options = optimset('Display','off');
switch P.case
    case 0 % do nothing
        P.label = 'true';
    case 1 % estimate a and b, unconstrained
        X   = [C 1+0*C];
        Y   = F;
        B   = X\Y;
        P.a = B(1);
        P.b = B(2);
        P.label = 'unconstrained';
    case 2 % estimate a, b > 0
        X   = [C 1+0*C];
        Y   = F;
        B   = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[inf inf],[P.a, P.b],options);
        P.a = B(1);
        P.b = B(2);
        P.label = 'a,b>0';
    case 3 % estimate a and b, 0<a<1, b>0
        X   = [C 1+0*C];
        Y   = F;
        B   = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[1 inf],[],options);
        P.a = B(1);
        P.b = B(2);
        P.label = '0<a<1, b>0';
    case 4 % estimate only a, unconstrained
        X   = C;
        Y   = F;
        P.a = X\Y;
        P.label = 'a unconstrained';
    case 5 % a>0
        X   = C;
        Y   = F;
        P.a = X\Y; if P.a<0, P.a=0; end
        P.label = 'a>0';
    case 6 % 0<a<1
        X   = C;
        Y   = F;
        P.a = X\Y; if P.a<0, P.a=0; elseif P.a>1, P.a=1; end
        P.label = '0<a<1';
    case 7 % estimate b > 0
        X   = F-P.a*C;
        Y   = 1+0*C;
        P.b = X\Y;
        if P.b<0, P.b=0; end
        P.label = 'b>0';
end

P.lam   = sum(n)/(T*dt);
DD      = norm(F-P.a*C-P.b)^2;
P.sig   = sqrt(DD/T);
l       = GetLik(F,n,T,dt,P);

end