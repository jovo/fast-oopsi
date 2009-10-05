function [n P]=FOOPSI2_54(F,P,Sim)
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
%
% Remarks on revisions:
% 1_7:  no longer need to define Sim.Plot (ie, if it is not defined, default
%       is to not plot, but if it is defined, one can either plot or not)
% 1_8:  cleaned up code from 1_7, and made Identity matrix outside loop,
%       which gets diagonal replaced inside loop (instead of calling speye in
%       loo)
% 1_9:  mean subtract and normalize max(F)=1 such that arbitrary scale and
%       offset shifts do no change results.
% 2:    removed normalize.  takes either a row or column vector.
%       doesn't require any Sim fields other than Sim.dt. also, we estimate
%       parameters now using FastParams code (which is the same as the one used
%       to estimate params given the real spikes, for debugging purposes)
% 2_1:  also estimate mu
% 2_2:  forgot to make this one :)
% 2_3:  fixed a bunch of bugs.  this version works to infer and learn, but
%       fixes mu in above model.
% 2_4:  to my knowledge, this one works, but requires fixing 'mu' and 'a' in
%       the above model. I also normalize between 0 and 1
% 2_41: reparameterized for stability.  uses constrained optimization. this
% works assuming nu=0 and rho=1.
% 2_42: works for arbitrary rho
% 2_43: fixed bugs so that M is only T-1 x T-1. cleaned up names and stuff.
% 2_431:made M TxT again
% 2_432:added baseline (in progress)
% 2_5:  dunno
% 2_51: removed rho and nu
% 2_52: a=alpha, b=beta, in code
% 2_53: threshold n s.t. n \in \{0,1\} before estimating parameters
% 2_54: allow for F to be a vector at each time step

%% initialize stuff
fprintf('\nFOOPSI2_54\n')

% get F "right"
siz     = size(F);                              % make sure F is a column vector
if siz(1)>1 && siz(2) >1
    error('F must be a vector')
elseif siz(1)==1 && siz(2)>1
    F=F';
end
Q=P;

% define some stuff for brevity
T       = Sim.T;                                % # of time steps
dt      = Sim.dt;                               % time step size
u       = 1/(2*P.sig^2);                        % scale of variance

% define some stuff for speed
O   = 1+0*F(:,1);                               % init a unity vector
M   = spdiags([-P.gam*O O], -1:0,T,T);          % matrix transforming calcium into spikes, ie n=M*C
I   = speye(T);                                 % create out here cuz it must be reused
H1  = 2*u*P.a^2*I;                                % pre-compute matrix for hessian
H2  = I;                                        % another one
d0  = 1:T+1:T^2;                                % index of diagonal elements of TxT matrices
d1  = 2:T+1:T^2;                                % index of off-diagonal elements (the diagonal below the diagonal) of TxT matrices

% if we are not estimating parameters
if ~isfield(Sim,'MaxIter') || Sim.MaxIter==0,
    Sim.MaxIter = 0;
    n = FastFilter(F,P);                        % infer approximate MAP spike train, given initial parameter estimates
    l = [];
else
    % initialize some stuff
    n   = O/P.lam;                              % spike train
    C   = filter(1,[1 -P.gam],n);               % calcium concentratin
    DD  = (F-C)'*(F-C);                         % squared error
    l   = zeros(Sim.MaxIter,1);                 % extize likelihood
    l(1)= .5*T*log(2*pi*P.sig^2) + DD/(2*P.sig^2) - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

    % prepare stuff for plotting
    if isfield(Sim,'Plot'), DoPlot = Sim.Plot; else DoPlot = 0; end
    if DoPlot == 1
        figure(104), clf
        %         fprintf('a=%.2f, b=%.2f, sig=%.2f, lam=%.2f\n',P.a, P.b, P.sig, P.gam, P.lam)
    end
end

%% here is the main loop
for i=1:Sim.MaxIter

    [n C]   = FastFilter(F,P);                  % infer approximate MAP spike train, given most recent parameter estimates
    P       = FastParams3_2(F,C,n,T,dt,Q);      % update parameters

%     if Sim.thresh==1
%         nsort   = sort(n);
%         nthr    = nsort(round(0.80*T));
%         nn      = n;
%         nn(nn<=nthr)=0;
%         nn(nn>nthr)=1;
%         CC = filter(1,[1 -P.gam],nn);
%         options = optimset('Display','off');
%         P   = Q;                                                % initialize all parameters
%         X   = [CC 1+0*C];
%         Y   = F;
%         B   = quadprog(X'*X,-X'*Y,[],[],[],[],[0 0],[inf inf],[Q.a, Q.b],options);
%         P.a = B(1);
%         P.b = B(2);
% 
%         DD      = norm(F-P.a*C-P.b)^2;
%         P.sig   = sqrt(DD/T);
%         P.lam   = sum(n)/(T*dt);
% 
%         c       = 1/(2*P.sig^2);                    % scale of variance
%         P.l     = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(nn);% initialize likelihood
%     else
%         P       = FastParams3_2(F,C,n,T,dt,Q);      % update parameters
%     end


    l(i+1) = P.l;                               % update likelihood and display stuff (if requested)
    if DoPlot == 1
        subplot(311), hold on, plot(i+1,l(i+1),'o'), axis('tight')
        subplot(312), cla, hold on, plot(F,'.k'), plot(P.a*C+P.b,'b'),  axis('tight')
        subplot(313), cla, bar(n/max(n),'EdgeColor','r','FaceColor','r'), axis('tight'), drawnow
        if Sim.thresh==1, hold on,
            plot(nthr/max(n)*(1+0*n)),
            hold off
        end
        fprintf('a=%.2f, b=%.2f, sig=%.2f, gam=%.2f, lam=%.2f, lik=%.2f\n',...
            P.a, P.b, P.sig, P.gam, P.lam, l(i+1))
    end

    if abs(l(i+1)-l(i))<1e-3, break, end        % stopping criterion
    %     if l(i+1)>l(i), display(l(i+1)-l(i)), keyboard, end
    if l(i+1)>l(i), display(l(i+1)-l(i)), end
end
%% main loop over
P.l = l(1:i+1);

    function [n C] = FastFilter(F,P)

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
%                 g   = -2*u*P.a*(F - P.a*C - P.b) + P.lam*dt*sumM - z*sum((M*C).^(-1)*ones(1,T).*M,1)';  % gradient
                g   = -2*u*P.a*(F - P.a*C - P.b) +P.lam*dt*sumM - z*M'*(n.^-1);  % gradient
                H2(d0) = n.^-2;
%                 U= ((M*C).^(-1)*ones(1,T)).* M;
%                 H   = H1 + z*U'*U;           % Hessian
                H   = H1 + z*(M'*H2*M);           % Hessian
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
    end

end