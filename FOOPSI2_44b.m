function [n P]=FOOPSI2_44(F,P,Sim)
% this function solves the following optimization problem:
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
%   gam:    "decay" (ie, tau=dt/(1-gam)
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
%
% Remarks on revisions:
% 1_7: no longer need to define Sim.Plot (ie, if it is not defined, default
% is to not plot, but if it is defined, one can either plot or not)
% 1_8: cleaned up code from 1_7, and made Identity matrix outside loop,
% which gets diagonal replaced inside loop (instead of calling speye in
% loo)
% 1_9: mean subtract and normalize max(F)=1 such that arbitrary scale and
% offset shifts do no change results.
%
% 2: removed normalize.  takes either a row or column vector.
% doesn't require any Sim fields other than Sim.dt. also, we estimate
% parameters now using FastParams code (which is the same as the one used
% to estimate params given the real spikes, for debugging purposes)
% 2_1: also estimate mu
% 2_2: forgot to make this one :)
% 2_3: fixed a bunch of bugs.  this version works to infer and learn, but
% fixes mu in above model.
% 2_4: to my knowledge, this one works, but requires fixing 'mu' and 'a' in
% the above model. I also normalize between 0 and 1
% 2_41: reparameterized for stability.  uses constrained optimization. this
% works assuming nu=0 and rho=1.
% 2_42: works for arbitrary rho
% 2_43: fixed bugs so that M is only T-1 x T-1
% 2_44: incorporate baseline
%%
fprintf('\nFOOPSI2_44\n')

% get F "right"
siz     = size(F);                              % make sure F is a column vector
if siz(1)>1 && siz(2) >1
    error('F must be a vector')
elseif siz(1)==1 && siz(2)>1                
    F=F';
end

% define some stuff for brevity
T       = length(F);                            % number of time steps
dt      = Sim.dt;                               % for brevity
c       = 1/(2*P.sig^2);                        % scale of variance

% original parameter estimate, keep these to provide as input into parameter estimation
Q.gam   = P.gam;
Q.nu    = P.nu;
Q.rho   = P.rho;

% define some stuff for speed
O   = 1+0*F(2:T);                               % init a unity vector
M   = spdiags([-P.gam*O O]/P.rho, 0:1,T-1,T-1); % matrix transforming calcium into spikes, ie n=M*C
M   = [-P.nu/P.rho*O M];                        % append vector for baseline
I   = speye(T);                                 % create out here cuz it must be reused
H1  = 2*c*I;                                    % pre-compute matrix for hessian
H2  = I(1:T-1,1:T-1);                           % another one
c1  = 1:T-1;                                    % 1st colomn of any TxT matrix (used for M)
d0  = 1:T:(T-1)^2;                              % index of diagonal elements
d1  = 500:T:(T-1)^2;                            % index of 1-off-diagonal elements
d2  = 2*T-1:T:(T-1)^2;                          % index of 2-off-diagonal elements

% if we are not estimating parameters
if ~isfield(Sim,'MaxIter') || Sim.MaxIter==0,
    Sim.MaxIter=0;
    [n C]   = FastFilter(F,P);                  % infer approximate MAP spike train, given initial parameter estimates
else
    % initialize some stuff
    n       = O/P.lam;                          % spike train
    C       = filter(P.rho,[1 -P.gam],n+P.nu/P.rho); % calcium concentratin
    DD      = (F(1:T-1)-C)'*(F(1:T-1)-C);       % squared error
    lik     = zeros(1,Sim.MaxIter);             % extize likelihood
    lik(1)  = .5*T*log(2*pi*P.sig^2) + DD/(2*P.sig^2) - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

    % prepare stuff for plotting
    if isfield(Sim,'Plot'), DoPlot = Sim.Plot; else DoPlot = 0; end
    if DoPlot == 1
        figure(104), clf
        fprintf('sig=%.2f, gam=%.2f, nu=%.2f, rho=%.2f, lam=%.2f, lik=%.2g\n',...
            P.sig, P.gam, P.nu, P.rho, P.lam, lik(1))
    end
end

for i=1:Sim.MaxIter

    [n C]   = FastFilter(F,P);                  % infer approximate MAP spike train, given most recent parameter estimates
    P       = FastParams2_7(F,C,n,T,dt,Q);      % update parameters

    lik(i+1) = P.lik;                           % update likelihood and display stuff (if requested)
    if DoPlot == 1
        subplot(311), hold on, plot(i+1,lik(i+1),'o'), axis('tight')
        subplot(312), cla, hold on, plot(F,'.k'), plot(C,'b'),  axis('tight')
        subplot(313), cla, bar(n,'EdgeColor','r','FaceColor','r'), axis('tight'), drawnow
        fprintf('sig=%.2f, gam=%.2f, nu=%.2f, rho=%.2f, lam=%.2f, lik=%.2g\n',...
            P.sig, P.gam, P.nu, P.rho, P.lam, lik(i+1))
    end

    if abs(lik(i+1)-lik(i))<1e-3, break, end    % stopping criterion
end

    function [n C] = FastFilter(F,P)

        e       = 1;                            % weight on barrier function
        c       = 1/(2*P.sig^2);                % scale of variance
        n       = O*(e/P.lam);                  % initialize spike train
        C       = [1; filter(P.rho,[1, -P.gam],n+P.nu/P.rho)]; % initialize calcium
        M(d1)   = -P.gam/P.rho;                 % update d1 of M
        M(d2)   = 1/P.rho;                      % matrix transforming calcium into spikes, ie n=M*C
        M(c1)   = -P.nu*P.rho;                  % update 1st col of M
        H1(d0)  = 2*c;                          % for expediency
        sumM    = sum(M)';                      % for expediency
        
        while e>1e-13                           % this is an arbitrary threshold

            D = F-C;                            % difference vector
            L = c*D'*D+P.lam*dt*sum(n)-e*sum(log(n));  % Likilihood function using C
            s = 1;                              % step size
            d = 1;                              % direction
            while norm(d)>5e-2 && s > 1e-3      % converge for this e (again, these thresholds are arbitrary)
                g   = -2*c*D + P.lam*dt*sumM - e*M'*(n.^-1);  % gradient
                H2(d0) = n.^-2;
                H   = H1 + 2*e*M'*H2*M;         % Hessian                
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
                    D   = F-C1;
                    L1  = c*D'*D+P.lam*dt*sum(n)-e*sum(log(n));
                    s   = s/2;                  % if step increases objective function, decrease step size
                end
                C = C1;                         % update C
                L = L1;                         % update L
            end
            e=e/10;                             % reduce e (sequence of e reductions is arbitrary)
        end
    end

end