function [n P]=FSOOPSI1_0(F,P,Sim)
% this function solves the following optimization problem:
% n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
%
% F_t = a*C_t + b + sigma*eps_t
% C_t = gam C_{t-1} + nu + rho*n_t
%
% where \{F_t,a,b,eps_t\} are vectors, each eps_{i,t} ~N(0,1), and n_t ~
% Exp(lam*dt).  we take an
% "interior-point" approach. each step with solved in O(T)
% time by utilizing gaussian elimination on the tridiagonal hessian, as
% opposed to the O(T^3) time typically required for non-negative
% deconvolution.
%
% Input----
% F:        fluorescence movie
% P.        structure of neuron parameters
%   alpha:  scale vector
%   beta:   offset vector
%   sig:    standard deviation
%   gam:  "decay" (ie, tau=dt/(1-gam)
%   nu:     "baseline" (ie, C_b = nu/(1-gam)
%   rho:    jump size
% Sim.      structure of simulation parameters
%   dt:     time step size
%   Plot:   whether to plot results        (if not set, default is no)
%   MaxIter:maximum number of iterations   (typically set to 50)
%
% Output---
% n:        inferred spike train
% P:        inferred parameter structure

%% initialize stuff
fprintf('\nFOOPSI2_53\n')
Q=P;

% define some stuff for brevity
T       = Sim.T;                            % number of time steps
dt      = Sim.dt;                               % time step size
u       = 1/(2*P.sig^2);                        % scale of variance

% define some stuff for speed
O   = 1+0*F(:,1);
M   = spdiags([-P.gam*O O], -1:0,T,T);% matrix transforming calcium into spikes, ie n=M*C
I   = speye(T);                                 % create out here cuz it must be reused
H1  = I;                                        % pre-compute matrix for hessian
H2  = I;                                        % another one
d0  = 1:T+1:T^2;                                % index of diagonal elements of TxT matrices
d1  = 2:T+1:T^2;                                % index of off-diagonal elements (the diagonal below the diagonal) of TxT matrices

[n C]   = FastFilter(F,P);                  % infer approximate MAP spike train, given most recent parameter estimates
P       = FastParams3_2(F,C,n,T,dt,Q);      % update parameters

    function [n C] = FastFilter(F,P)

        z       = 1;                            % weight on barrier function
        u       = 1/(2*P.sig^2);                % scale of variance
        n       = O*(z/P.lam);                  % initialize spike train
        C       = filter(1,[1, -P.gam],n);      % initialize calcium
        M(d1)   = -P.gam;                       % matrix transforming calcium into spikes, ie n=M*C
        sumM    = sum(M)';                      % for expediency
        H1(d0)  = 2*u*P.a*P.a';                    % for expediency

        while z>1e-13                           % this is an arbitrary threshold

            D = F-P.a'*C-P.b';                    % difference vector
            L = u*norm(D)^2+P.lam*dt*sum(n)-z*sum(log(n));  % Likilihood function using C
            L2 = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));  % Likilihood function using C
            s = 1;                              % step size
            d = 1;                              % direction
            a2= P.a'*P.a;                       
            while norm(d)>5e-2 && s > 1e-3      % converge for this z (again, these thresholds are arbitrary)
                g   = -2*u*P.a*(F - P.a*C - P.b) + P.lam*dt*sumM - z*M*(n.^-1);  % gradient
                g0   = -2*u*(F*P.a - a2*C-P.a*P.b*O)+ P.lam*dt*sumM - z*M*(n.^-1);  % gradient
                g1   = -2*u*P.a*(F - P.a*C - P.b) + P.lam*dt*sumM - z*M*(n.^-1);  % gradient
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
                    D   = F-P.a*C-P.b;
                    L1  = u*norm(D)^2+P.lam*dt*sum(n)-z*sum(log(n));
                    s   = s/2;                  % if step increases objective function, decrease step size
                end
                C = C1;                         % update C
                L = L1;                         % update L
            end
            z=z/10;                             % reduce z (sequence of z reductions is arbitrary)
        end
    end

end