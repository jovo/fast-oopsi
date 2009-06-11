function [n C] = FastFilter1_0(F,P)
% this function solves the following problem:
% 
% n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
%
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gam C_{t-1} + nu + rho*n_t, n_t ~ Poisson(n_t; p_t)
%
% and we approx the Poisson with an Exponential. we take an
% "interior-point" approach. each step with solved in O(T)
% time by utilizing gaussian elimination on the tridiagonal hessian, as
% opposed to the O(T^3) time typically required for non-negative
% deconvolution.
%
% Input----
% F:        fluorescence time series
% P.        structure of neuron parameters
%   alpha:  scale
%   beta:   offset
%   sig:    standard deviation
%   gam:  "decay" (ie, tau=dt/(1-gam)
%   nu:     "baseline" (ie, C_b = nu/(1-gam)
%   rho:    jump size
%
% Output---
% n:        inferred spike train
% P:        inferred parameter structure

z       = 1;                            % weight on barrier function
u       = 1/(2*P.sig^2);                % scale of variance
n       = O*(z/P.lam);                  % initialize spike train
C       = filter(1,[1, -P.gam],n);      % initialize calcium
M(d1)   = -P.gam;                       % matrix transforming calcium into spikes, ie n=M*C
sumM    = sum(M)';                      % for expediency
H1(d0)  = 2*u*P.a^2;                    % for expediency

while z>1e-13                           % this is an arbitrary threshold

    D = F-P.a*C-P.b;                    % difference vector
    L = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));  % Likilihood function using C
    s = 1;                              % step size
    d = 1;                              % direction
    while norm(d)>5e-2 && s > 1e-3      % converge for this z (again, these thresholds are arbitrary)
        g   = -2*u*P.a*(F - P.a*C - P.b) + P.lam*dt*sumM - z*M*(n.^-1);  % gradient
        H2(d0) = n.^-2;
        H   = H1 + z*M'*H2*M;           % Hessian
        d   = -H\g;                     % direction to step using newton-raphson
        hit = -n./(M*d);                % step within constraint boundaries
        hit(hit<0)=[];                  % ignore negative hits
        if any(hit<1)
            s = min(1,0.99*min(hit(hit>0)));
        else
            s = 1;
        end
        L1 = L+1;
        while L1>=L+1e-7                % make sure newton step doesn't increase objective
            C1  = C+s*d;
            n   = M*C1;
            D   = F-P.a*C1-P.b;
            L1  = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));
            s   = s/2;                  % if step increases objective function, decrease step size
        end
        C = C1;                         % update C
        L = L1;                         % update L
    end
    z=z/10;                             % reduce z (sequence of z reductions is arbitrary)
end