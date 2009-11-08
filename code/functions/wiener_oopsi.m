function [n P]=wiener_oopsi(F,dt,E)
% this function solves the following optimization problem:
% n = argmin sum_t ((F_t - C_t)^2 + lambda n_t^2
% where
% C_t = a C_{t-1} + n_t, where a=(1-dt/tau)
%
% and then estimates the parameters using:
%
% a = argmin_a ||Wa-Y||^2, where
%   W = -[C_1, C_{T-1}]
%   Y = [{-F_2 n_2}, {-F_T n_T}];
%
% the approach is called ridge regression.  it is essentially the least
% squares solution, regularized for sparsity.
%
% Input----
% F:    fluorescence time series
% dt:   time step size
% E.    structure of parameter estimates
%   tau:time constant
%   lam:prior weight = 1/(rate*A*dt)
%   sig:standard deviation of err
% 
% Output---
% n:    estimated spike train
% P:    estimated parameters

siz = size(F);
if siz(1)<siz(2), F=F'; end
F = F-mean(F);
F = F/max(abs(F));
T = length(F);                                  %# frames
o = 1+0*F;                                      %init a unity vector
M = spdiags([-E.gam*o o], -1:0,T,T);            %matrix transforming calcium into spikes, ie n=M*C
C = o;                                          %initialize calcium
P = E;
n = M*C;

DD1  = (F-C)'*(F-C);                            % squared error
DD2  = n'*n;                                    % squared error
lik  = .5*T*log(2*pi*P.sig^2) + DD1/(2*P.sig^2) + .5*T*log(2*pi*P.lam*dt) + DD2/(2*P.lam*dt); % initialize likelihood
old_lik = inf;
conv    = false;

while conv == false

    g = (C-F)/E.sig^2 + (((M*C)'*M)'-E.lam*dt*sum(M)')/(E.lam*dt);
    H = speye(T)/(E.sig^2) + M'*M/(E.lam*dt);
    C = C-H\g;
    N = M*C;

    old_lik = lik;
    lik = .5*T*log(2*pi*P.sig^2) + (F-C)'*(F-C)/(2*P.sig^2) + .5*T*log(2*pi*P.lam*dt) + (N'*N)/(2*P.lam*dt);
    if lik-old_lik<=-1e-4
        n = N;
        P = E;
        
        W = -C(1:end-1);
        Y = -F(2:end)+N(2:end);
        a = W\Y; %W'*Y/(W'*W);

        E.tau   = dt/(1-a);
        E.sig   = sqrt((F-C)'*(F-C)/T);
%         E.lam   = T*dt/sum(N);
    else
        conv = true;
    end
end