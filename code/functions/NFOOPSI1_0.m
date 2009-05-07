function [n P]=NFOOPSI1_0(F,P,Plin,Sim)
% this function solves the following optimization problem:
% n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
%
% F_t = alpha*S_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% S_t = C_t^n/(C_t^n+kd)
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
%   gam:    "decay" (ie, tau=dt/(1-gam)
%   nu:     "baseline" (ie, C_b = nu/(1-gam)
%   rho:    jump size
% Plin.     structure of neuron parameters for the linear case
% Sim.      structure of simulation parameters
%   dt:     time step size
%   Plot:   whether to plot results        (if not set, default is no)
%   MaxIter:maximum number of iterations   (typically set to 50)
%
% Output---
% n:        inferred spike train
% P:        inferred parameter structure

%% initialize stuff
fprintf('\nNFOOPSI1_0\n')

% get F "right"
siz     = size(F);                              % make sure F is a column vector
if siz(1)>1 && siz(2) >1
    error('F must be a vector')
elseif siz(1)==1 && siz(2)>1
    F=F';
end

% define some stuff for brevity
T   = Sim.T;                                    % # of time steps
dt  = Sim.dt;                                   % time step size
O   = 1+0*F(:,1);                               % init a unity vector
I   = speye(T);                                 % create out here cuz it must be reused
H1  = I;                                        % initialize memory for Hessian matrix
H2  = I;                                        % another one
d0  = 1:T+1:T^2;                                % index of diagonal elements of TxT matrices
d1  = 2:T+1:T^2;                                % index of off-diagonal elements (the diagonal below the diagonal) of TxT matrices

% define some stuff for speed
M   = spdiags([-P.gam*O O]/P.rho, -1:0,T,T);    % matrix transforming calcium into spikes, ie n=M*C
ML  = spdiags([-Plin.gam*O O], -1:0,T,T);       % matrix transforming calcium into spikes, ie n=M*C

n = NonlinFastFilter(F,P);                      % infer approximate MAP spike train, given initial parameter estimates

    function [n C] = NonlinFastFilter(F,P)

        z       = 1;                            % weight on barrier function
        u       = 1/(2*P.sig^2);                % scale of variance
        n       = O*(z/P.lam);                  % initialize spike train
        C       = filter(1,[1, -P.gam],n*P.rho);% initialize calcium
        S       = C./(C+P.k_d);
        M(d0)   = 1/P.rho;
        M(d1)   = -P.gam/P.rho;                 % matrix transforming calcium into spikes, ie n=M*C
        sumM    = sum(M)';                      % for expediency
        H1(d0)  = 2*u*P.a^2;                    % for expediency

%         uL      = 1/(2*Plin.sig^2);             % scale of variance
%         CL      = filter(1,[1, -Plin.gam],n);   % initialize calcium
%         ML(d0)  = 1;
%         ML(d1)  = -Plin.gam;                    % matrix transforming calcium into spikes, ie n=M*C
%         sumML   = sum(ML)';                     % for expediency
%         H1L(d0) = 2*u*Plin.a^2;                 % for expediency

        while z>1e-13                           % this is an arbitrary threshold

%             D = F-P.a*C-P.b;                    % difference vector
%             L = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));  % Likilihood function using C
%             s = 1;                              % step size
%             d = 1;                              % direction
%             while norm(d)>5e-2 && s > 1e-3      % converge for this z (again, these thresholds are arbitrary)
%                 g   = -2*u*P.a*(F - P.a*C - P.b) +P.lam*dt*sumM - z*M'*(n.^-1);  % gradient
%                 H2(d0) = n.^-2;
%                 H   = H1 + z*(M'*H2*M);         % Hessian
%                 d   = -H\g;                     % direction to step using newton-raphson
%                 hit = -n./(M*d);                % step within constraint boundaries
%                 hit(hit<0)=[];                  % ignore negative hits
%                 if any(hit<1)
%                     s = min(1,0.99*min(hit(hit>0)));
%                 else
%                     s = 1;
%                 end
%                 L1 = L+1;
%                 while L1>=L+1e-7                % make sure newton step doesn't increase objective
%                     C1  = C+s*d;
%                     n   = M*C1;
%                     D   = F-P.a*C1-P.b;
%                     L1  = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));
%                     s   = s/2;                  % if step increases objective function, decrease step size
%                 end
%                 C = C1;                         % update C
%                 L = L1;                         % update L
%             end
         
%             S = C./(C+P.k_d);           
            D = F-P.a*S-P.b;                    % difference vector
            L = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));  % Likilihood function using C
            s = 1;                              % step size
            d = 1;                              % direction
            while norm(d)>5e-2 && s > 1e-3      % converge for this z (again, these thresholds are arbitrary)
                D   = F-P.a*S-P.b;
                g   =  -2*D*P.a*P.k_d.*(C+P.k_d).^-2 ...
                    +P.lam*dt*sumM - z*M'*(n.^-1);  % gradient
                H2(d0) = n.^-2;
                x   =(-P.a*P.k_d-2*(C+P.k_d).*D).*((C+P.k_d).^-4);
                H1(d0) = x;
                H   = H1 + z*(M'*H2*M);         % Hessian
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
                    S   = C1./(C1+P.k_d);
                    D   = F-P.a*S-P.b;
                    L1  = u*D'*D+P.lam*dt*sum(n)-z*sum(log(n));
                    s   = s/2;                  % if step increases objective function, decrease step size
                end
                C = C1;                         % update C
                L = L1;                         % update L
            end
            z=z/10;                             % reduce z (sequence of z reductions is arbitrary)
        end
    end

end