function [n_best P_best]=FOOPSI2(F,P,Sim)
% this function solves the following optimization problem:
% n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
%
% n_t ~ Poisson(n_t; p_t)
% C_t = a C_{t-1} + n_t, where a=(1-dt/tau), and A is fixed at 1.
% F_t = C_t + eps_t, eps_t ~ N(mu,sig^2)
%
% and we approx the Poisson with an Exponential yielding
% n_best = argmin_{n: n_t >= 0 forall t} sum_t ((F_t - C_t)^2 + lambda n_t
%
% instead of solving this directly, we iteratively solve
%
% n_eta = argmin sum_t ((F_t - C_t)^2 + lambda n_t -eta*log(n_t)
%
% starting with eta at some reasonably large value, and then decreasing.
% note that -log(n_t) is a "barrier term" which softens the constraints
% mapping the problem into a concave one.  each step with solved in O(T)
% time by utilizing gaussian elimination on the tridiagonal hessian, as
% opposed to the O(T^3) time typically required for non-negative
% deconvolution.
%
% Input----
% F:    fluorescence time series
% P.    structure of neuron parameters
%   tau: time constant size
%   lam: prior weight = 1/(rate*dt)
%   sig: standard deviation of err
%   mu : mean of err
% Sim.  structure of simulation parameters
%   dt      : time step size
%   Plot    : whether to plot results        (if not set, default is no)
%   MaxIter : maximum number of iterations   (typically set to 50)
% Output---
% n_best:    inferred spike train
% P_best:    inferred parameter structure
%
% Remarks on revisions:
% 1_7: no longer need to define Sim.Plot (ie, if it is not defined, default
% is to not plot, but if it is defined, one can either plot or not)
%
% 1_8: cleaned up code from 1_7, and made Identity matrix outside loop,
% which gets diagonal replaced inside loop (instead of calling speye in
% loo)
%
% 1_9: mean subtract and normalize max(F)=1 such that arbitrary scale and
% offset shifts do no change results.
%
% 2_0: removed normalize.  takes either a row or column vector.
% doesn't require any Sim fields other than Sim.dt. also, we estimate
% parameters now using FastParam1 code (which is the same as the one used
% to estimate params given the real spikes, for debugging purposes)
% 
% 2_1: also estimate mu


%%
fprintf('\nFOOPSI2\n')

% make sure F is a column vector
siz     = size(F);                          
if siz(1)>1 && siz(2) >1
    error('F must be a vector')
elseif siz(1)==1 && siz(2)>1                
    F=F';
end

% determine number of interations, if not declared, set to 0
if ~isfield(Sim,'MaxIter'), Sim.MaxIter=0; end 

% define some stuff for brevity
T       = length(F);                        % number of time steps
dt      = Sim.dt;                           % for brevity
a       = 1 - dt/P.tau;                     % decay factor
c       = 1/(2*P.sig^2);                    % scale of variance
mu      = P.mu;                             % brevity

% define some stuff for speed
o       = 1+0*F;                            % init a unity vector
M       = spdiags([-a*o o], -1:0,T,T);      % matrix transforming calcium into spikes, ie n=M*C
I       = speye(T);                         % create out here cuz it must be reused
Hmat1   = 2*c*I;                            % pre-compute matrix for hessian
Hmat2   = I;                                % another one
diags   = 1:T+1:T^2;

% initialize some stuff
n       = o/P.lam;                          % spike train
C       = filter(1,[1 -(1-dt/P.tau)],n);    % calcium concentratin
DD      = (F-C-mu)'*(F-C-mu);                     % squared error
lik     = zeros(1,Sim.MaxIter);             % extize likelihood
lik(1)  = .5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood
minlik  = lik(1);                           % minimum likelihood achived so far

% prepare plotting stuff
DoPlot  = isfield(Sim,'Plot');
if DoPlot == 1 && Sim.Plot== 1
    figure(104), clf
    %     subplot(311), hold on, plot(1,lik(1),'o'), axis('tight')
    fprintf('lam=%.2f, tau=%.2f, sig=%.2f, lik=%.2f\n',P.lam,P.tau,P.sig,lik(1))
end

% if not estimating parameters
[n C]   = FastFilter(F);
n_best  = n;
P_best  = P;

% else iterate until achieving the max number of iterations (or lik stops
% changing)
for i=1:Sim.MaxIter

    [n C]   = FastFilter(F);
    P       = FastParams1(F,C,n,T,dt);

    %     % update tau
    %     W = C(1:end-1);
    %     Y = F(2:end)-n(2:end);
    %     a = W\Y;
    %     if a>1, a=1; elseif a<0, a=0; end           % make sure 'a' is within bounds
    %     M     = spdiags([-a*o o], -1:0,T,T);        % matrix transforming calcium into spikes, ie n=M*C
    %     P.tau = dt/(1-a);
    %
    %     % update sig
    %     DD      = (F-C)'*(F-C);                     % squared error
    %     P.sig   = sqrt(DD/T);
    %     c       = 1/(2*P.sig^2);
    %     Hmat1  = 2*c*I;
    %
    %     % update lambda
    %     P.lam = sum(n)/(T*dt);
    %
    %     lik(i+1) = .5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);

    lik(i+1) = P.lik;
    if DoPlot == 1 && Sim.Plot== 1
        subplot(311), hold on, plot(i+1,lik(i+1),'o'), axis('tight')
        subplot(312), cla, hold on, plot(F,'.k'), plot(C,'b'),  axis('tight')
        subplot(313), cla, bar(n,'EdgeColor','r','FaceColor','r'), axis('tight'), drawnow
        fprintf('lam=%.2f, tau=%.2f, sig=%.2f, lik=%.2f\n',P.lam,P.tau,P.sig,lik(i+1))
    end

    % the next line of code is only necessary because likelihood
    % doesn't increase with each step
    if lik(i+1)<minlik,
        minlik=lik(i+1); P_best = P; n_best = n;
    end

    % stopping criterion
    if abs(lik(i+1)-lik(i))<1e-3, break, end
end
P_best.i = i;


    function [n C] = FastFilter(F)

        eta = 1;                                % weight on barrier function
        n=o*(eta/P.lam);
        C=filter(1,[1, -a],n);

        while eta>1e-13                         % this is an arbitrary threshold

            D = F-C;                            % difference vector
            L = c*D'*D+P.lam*dt*sum(n)-eta*sum(log(n));  % Likilihood function using C
            s = 1;                              % step size
            d = 1;                              % direction
            while norm(d)>5e-2 && s > 1e-3      % converge for this eta (again, these thresholds are arbitrary)
                g   = -2*c*D + P.lam*dt*sum(M)' - eta*M'*(n.^-1);  % gradient
                Hmat2(diags) = n.^-2;
                H   = Hmat1 + 2*eta*M'*Hmat2*M; % Hessian
                d   = -H\g;                     % direction to step using newton-raphson
                hit = -n./(M*d);                % step within constraint boundaries
                hit(hit<0)=[];                  % ignore negative hits
                if any(hit<1)
                    s = min(1,0.99*min(hit(hit>0)));
                else
                    s = 1;
                end
                L_new = L+1;
                while L_new>=L+1e-7             % make sure newton step doesn't increase objective
                    C_new   = C+s*d;
                    n       = M*C_new;
                    D       = F-C_new;
                    L_new   = c*D'*D+P.lam*dt*sum(n)-eta*sum(log(n));
                    s       = s/2;              % if step increases objective function, decrease step size
                end

                C = C_new;                      % update C
                L = L_new;                      % update L

            end
            eta=eta/10;                         % reduce eta (sequence of eta reductions is arbitrary)
        end
    end

end