function [n P]=FOOPSI2_42(F,P,Sim)
% this function solves the following optimization problem:
% n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
%
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma C_{t-1} + nu + rho*n_t, n_t ~ Poisson(n_t; p_t)
%
% and we approx the Poisson with an Exponential. we take an
% "interior-point" approach. each step with solved in O(T)
% time by utilizing gaussian elimination on the tridiagonal hessian, as
% opposed to the O(T^3) time typically required for non-negative
% deconvolution.
%
% Input----
% F:    fluorescence time series
% P.    structure of neuron parameters
%   alpha:  scale
%   beta:   offset
%   sig:    standard deviation
%   gamma:  "decay" (ie, tau=dt/(1-gamma)
%   nu:     "baseline" (ie, C_b = nu/(1-gamma)
%   rho:    jump size
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
% 2: removed normalize.  takes either a row or column vector.
% doesn't require any Sim fields other than Sim.dt. also, we estimate
% parameters now using FastParams code (which is the same as the one used
% to estimate params given the real spikes, for debugging purposes)
%
% 2_1: also estimate mu
% 2_2: forgot to make this one :)
% 2_3: fixed a bunch of bugs.  this version works to infer and learn, but
% fixes mu in above model.
% 2_4: to my knowledge, this one works, but requires fixing 'mu' and 'a' in
% the above model. I also normalize between 0 and 1
% 2_41: reparameterized for stability.  uses constrained optimization. this
% works assuming nu=0 and rho=1.
% 2_42: works for arbitrary rho
%%
fprintf('\nFOOPSI2_42\n')

% get F "right"
% F       = F - min(F);                       % min subtraction
% F       = F / max(F);                       % normalize so max(F)=1
siz     = size(F);                          % make sure it is a vector
if siz(1)>1 && siz(2) >1
    error('F must be a vector')
elseif siz(1)==1 && siz(2)>1                % make sure it is a column vector
    F=F';
end

% define some stuff for brevity
T       = length(F);                        % number of time steps
dt      = Sim.dt;                           % for brevity
c       = 1/(2*P.sig^2);                    % scale of variance

% original parameter estimate, keep these to provide as input into
% parameter estimation
Q.gamma = P.gamma;
Q.nu    = P.nu;
Q.rho   = P.rho;

% define some stuff for speed
O       = 1+0*F;                            % init a unity vector
M       = spdiags([-P.gamma*O O]/P.rho, -1:0,T,T); % matrix transforming calcium into spikes, ie n=M*C
I       = speye(T);                         % create out here cuz it must be reused
Hmat1   = 2*c*I;                            % pre-compute matrix for hessian
Hmat2   = I;                                % another one
diags   = 1:T+1:T^2;                        % index of diagonal elements of TxT matrices
offdiags=2:T+1:T^2;                         % index of off-diagonal elements (the diagonal below the diagonal) of TxT matrices

% if we are not estimating parameters
if ~isfield(Sim,'MaxIter') || Sim.MaxIter==0,
    Sim.MaxIter=0;
    [n C]   = FastFilter(F,P);
else
    % initialize some stuff
    n       = O/P.lam;                          % spike train
    C       = filter(P.rho,[1 -P.gamma],n);    % calcium concentratin
    DD      = (F-C)'*(F-C);                     % squared error
    lik     = zeros(1,Sim.MaxIter);             % extize likelihood
    lik(1)  = .5*T*log(2*pi*P.sig^2) + DD/(2*P.sig^2) - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

    % prepare stuff for plotting
    if isfield(Sim,'Plot'), DoPlot = Sim.Plot; else DoPlot = 0; end
    if DoPlot == 1
        figure(104), clf
        fprintf('sig=%.2f, gamma=%.2f, nu=%.2f, rho=%.2f, lam=%.2f, lik=%.2g\n',P.sig, P.gamma, P.nu, P.rho, P.lam, lik(1))
    end
end

for i=1:Sim.MaxIter

    [n C]   = FastFilter(F,P);
    P       = FastParams2_6(F,C,n,T,dt,Q);

    lik(i+1) = P.lik;
    if DoPlot == 1
        subplot(311), hold on, plot(i+1,lik(i+1),'o'), axis('tight')
        subplot(312), cla, hold on, plot(F,'.k'), plot(C,'b'),  axis('tight')
        subplot(313), cla, bar(n,'EdgeColor','r','FaceColor','r'), axis('tight'), drawnow
        fprintf('sig=%.2f, gamma=%.2f, nu=%.2f, rho=%.2f, lam=%.2f, lik=%.2g\n',P.sig, P.gamma, P.nu, P.rho, P.lam, lik(i+1))
    end

    % stopping criterion
    if abs(lik(i+1)-lik(i))<1e-3, break, end
end

    function [n C] = FastFilter(F,P)

        e = 1;                                % weight on barrier function
        c   = 1/(2*P.sig^2);                    % scale of variance
        n   = O*(e/P.lam);                    % initialize spike train
        C   = filter(P.rho,[1, -P.gamma],n);              % initialize calcium
        M(offdiags) = -P.gamma/P.rho;                       % matrix transforming calcium into spikes, ie n=M*C
        M(diags)    = 1/P.rho;
        Hmat1(diags)= 2*c;                      % pre-compute matrix for hessian
        sumM        = sum(M)';
        
        while e>1e-13                         % this is an arbitrary threshold

            D = F-C;                         % difference vector
            L = c*D'*D+P.lam*dt*sum(n)-e*sum(log(n));  % Likilihood function using C
            s = 1;                              % step size
            d = 1;                              % direction
            while norm(d)>5e-2 && s > 1e-3      % converge for this e (again, these thresholds are arbitrary)
                g   = -2*c*D + P.lam*dt*sumM - e*M'*(n.^-1);  % gradient
                Hmat2(diags) = n.^-2;
                H   = Hmat1 + 2*e*M'*Hmat2*M; % Hessian                
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
                    L_new   = c*D'*D+P.lam*dt*sum(n)-e*sum(log(n));
                    s       = s/2;              % if step increases objective function, decrease step size
                end
                C = C_new;                      % update C
                L = L_new;                      % update L
            end
            e=e/10;                         % reduce e (sequence of e reductions is arbitrary)
        end
    end

end