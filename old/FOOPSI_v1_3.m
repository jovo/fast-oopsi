function [n_best P_best]=FOOPSI_v1_3(F,P,Sim)
% this function solves the following optimization problem:
% n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
% n_t ~ Binomial(n_t; p_t)
% C_t = a C_{t-1} + A n_t, where a=(1-dt/tau), and A is fixed at 1.
% F_t = C_t + eps_t, eps_t ~ N(0,sig^2)
% and a little simplification, we have
% n_best = argmin sum_t ((F_t - C_t)^2 + lambda n_t
% with constraints: n_t>0 forall t, and
%
% instead of solving this directly, we iteratively solve
%
% n_eta = argmin sum_t ((F_t - C_t)^2 + lambda n_t -eta*log(n_t)
%
% starting with eta at some reasonably large value, and then decreasing.
% note that -log(n_t) is a "barrier term" which softens the constraints
% mapping the problem into a concave one.  each step with solved in O(T)
% time by utilizing gaussian elimination on the tridiagonal hessian.
%
% Input----
% F:    fluorescence time series
% P.    structure of neuron parameters
%   tau: time constant size
%   lam: prior weight = 1/(rate*A*dt)
%   sig: standard deviation of err
% Sim.  structure of simulation parameters
%   dt      : time step size
%   T       : number of time steps
%   Mstep   : whether to estimate parameters (typically set to yes)
%   Plot    : whether to plot results        (typically set to no)
%   MaxIter : maximum number of iterations   (typically set to 50)
% Output---
% n_best:    inferred spike train
% P_best:    inferred parameter structure
%
% note that in this version we assume that we know lambda, so we do not
% estimate it. it also assumes obsevations at every time step. and it
% assumes a 1D input corresponding to the spatially averaged fluorescence
% signal.

clc, fprintf('\nFOOPSI_v1_3\n')
T       = Sim.T;                            % # of time steps
a       = 1 - Sim.dt/P.tau;                 % decay factor
o       = 1+0*F;                            % init a unity vector
M       = spdiags([-a*o o], -1:0,T,T);      % matrix transforming calcium into spikes, ie n=M*C
Tmat    = speye(T);                         % create out here cuz it must be reused
c       = 1/(2*P.sig^2);                    % scale of variance
Hmat    = 2*c*Tmat;                         % Hessian
eta     = 1;                                % weight on barrier function
n       = o*eta/P.lam;                      % spike train
C       = filter(1,[1 -(1-Sim.dt/P.tau)],n);% calcium concentratin
lik     = zeros(1,Sim.MaxIter);             % extize likelihood
lik(1)  = c*(F-C)'*(F-C)+P.lam*sum(n)+eta*sum(log(n));% initialize likelihood
i       = 2;                                % iteration number
minlik  = inf;                              % minimum likelihood achived so far
conv    = false;                            % param estimate converge?
figure(100), clf
while conv == false

    eta = 1;                                % weight on barrier function
    n=o*(eta/P.lam);
    C=filter(1,[1 -a],n);

    while eta>1e-13                         % this is an arbitrary threshold

        D = F-C;                            % difference vector
        n = M*C;                            % get spike train assuming a particular C
        L = c*D'*D+P.lam*sum(n)-eta*sum(log(n));  % Likilihood function using C
        s = 1;                              % step size
        d = 1;                              % direction
        while norm(d)>5e-2 && s > 1e-3     % converge for this eta (again, these thresholds are arbitrary)
            g   = -2*c*D+P.lam*Sim.dt*sum(M)'-eta*M'*(n.^-1);  % gradient
            H   = Hmat+2*eta*M'*spdiags(n.^-2,0,T,T)*M; % Hessian
            d   = -H\g;                     % direction to step using newton-raphson
            hit = -n./(M*d);                % step within constraint boundaries
            hit(hit<0)=[];                  % ignore negative hits
            if any(hit<1)
                s = min(1,0.99*min(hit(hit>0)));
            else
                s = 1;
            end
            L_new = L+1;
            while L_new>=L+1e-2             % make sure newton step doesn't increase objective
                C_new   = C+s*d;
                n       = M*C_new;
                D       = F-C_new;
                L_new   = c*D'*D+P.lam*sum(n)-eta*sum(log(n));
                s       = s/2;              % if step increases objective function, decrease step size
            end

            C = C_new;                      % update C
            L = L_new;                      % update L

        end
        eta=eta/10;                         % reduce eta (sequence of eta reductions is arbitrary)
    end

    if Sim.Mstep==1;
        lik(i)     = L;
        if i<Sim.MaxIter % lik(i)-lik(i-1)<=-1e-4
            if Sim.Plot == true
                subplot(311), hold on, plot(i,lik(i),'o'), axis('tight')
                subplot(312), cla, hold on, plot(F,'.k'), plot(C,'b'),  axis('tight')
                subplot(313), cla, bar(n,'EdgeColor','r','FaceColor','r'), axis('tight'), drawnow
                fprintf('tau=%.2f, sig=%.2f, lam=%.2f, L=%.2f\n',P.tau,P.sig,P.lam,lik(i-1))
            end

            % the next line of code is only necessary because likelihood
            % doesn't increase with each step
            if lik(i)<minlik, minlik=lik(i); P_best = P; n_best = n; end

            i = i+1;

            % update tau
            %             W = C(1:end-1);
            %             Y = F(2:end)-n(2:end);
            %             a = W'*Y/(W'*W);
            %             if a>1, a=1; elseif a<0, a=0; end           % make sure 'a' is within bounds
            %             M     = spdiags([-a*o o], -1:0,T,T);        % matrix transforming calcium into spikes, ie n=M*C
            %             P.tau = Sim.dt/(1-a);

            % update sig
            P.sig = sqrt((F-C)'*(F-C)/T);
            c     = 1/(2*P.sig^2);
            Hmat  = 2*c*Tmat;

            % update lambda
            %             P.lam = T/(Sim.dt*sum(n));
        else
            conv = true;
        end
    else
        n_best = n;
        P_best = P;
        conv = true;
    end
end