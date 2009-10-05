function [P l] = FastParams4_0(F,C,n,T,dt,P)
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
X       = [C 1+Z(1:T)];
for ii=1:Sim.Np
    Y   = F(:,ii);
    B   = X\Y;
    for j=1:Nc
        P.a(ii,j) = B(j);
    end
    P.b(ii) = B(end);
end
nnorm   = n./max(n);
P.lam   = sum(nnorm)'/(T*dt);
P.sig   = sqrt(DD/T);

end