function [n_best P_best]=FOOPSI_v3_03_01(F,P,Meta,User)
% this function solves the following optimization problem:
% (*) n_best = argmax_{n >= 0} P(n | F)
% which is a MAP estimate for the most likely spike train given the
% fluorescence signal.  given the model:
%
% C_t = gam C_{t-1} + nu + rho*n_t, n_t ~ Poisson(n_t; p_t)
% F_t = alpha*C_t + beta + sigma*eps_t, eps_t ~ N(0,1)
%
% if F_t is a vector, then alpha and beta are BOTH vectors as well
% we approx the Poisson with an Exponential. we take an
% "interior-point" approach to impose the nonnegative contraint on (*). each step with solved in O(T)
% time by utilizing gaussian elimination on the tridiagonal hessian, as
% opposed to the O(T^3) time typically required for non-negative
% deconvolution.
%
% Input----
% F:        fluorescence time series (can be a vector (1 x T) or a matrix (Np x T)
% P.        structure of neuron parameters
%   alpha:  scale
%   beta:   offset
%   sig:    standard deviation
%   gam:    decayish (ie, tau=dt/(1-gam)
%   lam:    firing rate-ish
% Meta.      structure of simulation parameters
%   dt:     time step size
%   T:      # of time steps
%   h:      height of ROI (assumes square ROI) (# of pixels)
%   w:      width of ROI (assumes square ROI) (# of pixels)
%   Np:     # of pixels in ROI
% User.     structure of User defined parameters
%   Plot:   whether to plot results        (if not set, default is no)
%   MaxIter:maximum number of iterations   (typically set to 50)
%   Nc:     # of cells within ROI
%   Thresh: whether to threshold infered spike train before updating parameters
%   n:      if true spike train is known, and we are plotting, plot it

% Output---
% n:        inferred spike train
% P:        inferred parameter structure
%
% Remarks on revisions:
% 1_7:      no longer need to define Meta.Plot (ie, if it is not defined, default
%           is to not plot, but if it is defined, one can either plot or not)
% 1_8:      cleaned up code from 1_7, and made Identity matrix outside loop,
%           which gets diagonal replaced inside loop (instead of calling speye in
%           loo)
% 1_9:      mean subtract and normalize max(F)=1 such that arbitrary scale and
%           offset shifts do no change results.
% 2:        removed normalize.  takes either a row or column vector.
%           doesn't require any Meta fields other than Meta.dt. also, we estimate
%           parameters now using FastParams code (which is the same as the one used
%           to estimate params given the real spikes, for debugging purposes)
% 2_1:      also estimate mu
% 2_2:      forgot to make this one :)
% 2_3:      fixed a bunch of bugs.  this version works to infer and learn, but
%           fixes mu in above model.
% 2_4:      to my knowledge, this one works, but requires fixing 'mu' and 'a' in
%           the above model. I also normalize between 0 and 1
% 2_41:     reparameterized for stability.  uses constrained optimization. this
%           works assuming nu=0 and rho=1.
% 2_42:     works for arbitrary rho
% 2_43:     fixed bugs so that M is only T-1 x T-1. cleaned up names and stuff.
% 2_431:    made M TxT again
% 2_432:    added baseline (in progress)
% 2_5:      dunno
% 2_51:     removed rho and nu
% 2_52:     a=alpha, b=beta, in code
% 2_53:     threshold n s.t. n \in \{0,1\} before estimating parameters
% 2_54:     allow for F to be a vector at each time step
% 2_55:     fixed bugs, back to scalar F_t
% 2_56:     back to vector case, but no param estimate
% 2_57:     estimate spatial filter as well
% 2_58:     multiple cells (buggy)
% 2_59:     multiple cells, estimate {a,b}
% 3_01_01:  cleaning up a bit
% 3_02_01:  added input structure 'U' to control parameters that are 'User defined'
% 3_02_02:  don't need to include U in input, default values are set, and
%           rearranged some code, added some comments
% 3_02_03:  no more GetLik function, just inline, also plot true n if
%           available from User structure, and plot max lik
% 3_03_01:  made background a scalar

%% initialize stuff

% default user defined variables
U.Plot=1;
U.MaxIter=25;
U.Nc=1;
U.Thresh=1;

% if none or only some are defined by user, use defaults
if nargin == 3
    User=U;
else
    if isfield(User,'Plot'), else User.Plot=U.Plot; end
    if isfield(User,'MaxIter'), else User.MaxIter=U.MaxIter; end
    if isfield(User,'Nc'), else User.Nc=U.Nc; end
    if isfield(User,'Thresh'), else User.Thresh=U.Thresh; end
end

% define some stuff for brevity
Nc      = User.Nc;                              % # of cells
T       = Meta.T;                               % # of time steps
dt      = Meta.dt;                              % time step size
e       = 1/(2*P.sig^2);                        % scale of variance

% define some stuff for speed
Z   = zeros(Nc*T,1);                            % zero vector
M   = spdiags([repmat(-P.gam,T,1) repmat(Z,1,Nc-1) (1+Z)], -Nc:0,Nc*T,Nc*T);  % matrix transforming calcium into spikes, ie n=M*C
I   = speye(Nc*T);                              % create out here cuz it must be reused
H1  = I;                                        % initialize memory for Hessian matrix
H2  = I;                                        % another one
d0  = 1:Nc*T+1:(Nc*T)^2;                        % index of diagonal elements of TxT matrices
d1  = 1+Nc:Nc*T+1:(Nc*T)*(Nc*(T-1));            % index of diagonal elements of TxT matrices
l   = Z(1:User.MaxIter);                        % initialize likelihood
n   = FastFilter(F,P);                      % infer approximate MAP spike train, given initial parameter estimates
l(1)    = -inf;
l_max   = l(1);                                 % maximum likelihood achieved so far
n_best  = n;                                    % best spike train
P_best  = P;                                    % best parameter estimate

for i=2:User.MaxIter
    
    % update inferred spike train
    [n C DD]   = FastFilter(F,P);

    if min(n(:))<0, 
        disp('some shit is fucked'), keyboard
    end
    
    % generate regressor for spatial filter
    if User.Thresh==1
        CC=0*C;
        for j=1:Nc
            nsort   = sort(n(:,j));
            nthr    = nsort(round(0.95*T));
            nn      = Z(1:T);
            nn(n(:,j)<=nthr)=0;
            nn(n(:,j)>nthr)=1;
            CC(:,j) = filter(1,[1 -P.gam(j)],nn) + (1-P.gam(j))*P.b(j);
        end
    else
        CC      = C; 
    end

    % update spatial filter
    for ii=1:Meta.Np
        Y   = F(:,ii);
        B   = CC\Y;
        for j=1:Nc
            P.a(ii,j) = B(j);
        end
    end

    % estimate other parameters
    P.sig   = sqrt(DD/T);
    nnorm   = n./repmat(max(n),Meta.T,1);
    P.lam   = sum(nnorm)'/(T*dt);
        
    % update likelihood and keep results if they improved
    lik     = -Meta.T*Meta.Np*log(2*pi*P.sig^2)/2 -1/(2*P.sig^2)*DD;
    prior   = Meta.T*sum(P.lam*Meta.dt) - Meta.dt*P.lam'*sum(n)';
    l(i)    = lik + prior;

     % if this is the best one, keep n and P
    if l(i)>l_max                              
        n_best  = n;
        P_best  = P;
        l_max   = l(i);
    end

    % if lik doesn't change much (relatively), stop iterating
    if abs((l(i)-l(i-1))/l(i))<1e-5 || l(i-1)-l(i)>1e5; break; end 
    
    % plot results from this iteration
    if User.Plot == 1
        figure(400), nrows=Nc;                % plot spatial filter
        for j=1:Nc, subplot(1,nrows,j),
            imagesc(reshape(z1(P.a(:,j)),Meta.h,Meta.w)),
            title('a')
        end

        figure(401), clf, ncols=Nc+1;
        for j=1:Nc                              % plot inferred spike train
            h(j)=subplot(ncols,1,j); 
            bar(z1(n(:,j)))
            
            if isfield(User,'n'), hold on,
                stem(User.n(:,j),'LineStyle','none','Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',2)
            end
            title(['iteration ' num2str(i)]),
            set(gca,'XTickLabel',[])
            axis('tight')
        end
        subplot(ncols,1,ncols), plot(l(2:i))    % plot record of likelihoods
        title(['max lik ' num2str(l_max,4), ',   lik ' num2str(l(i),4)])
        set(gca,'XTickLabel',[])
        linkaxes(h,'xy')
        drawnow
    end
    
    % play sound to indicate iteration is over
    sound(3*sin(linspace(0,90*pi,2000)))        

end

P_best.l=l(1:i);                                % keep record of likelihoods for record

    function [n C DD] = FastFilter(F,P)

        % initialize n and C
        z = 1;                                  % weight on barrier function
        e = 1/(2*P.sig^2);                      % scale of variance
        n = repmat(z./P.lam,T,1);               % initialize spike train
        C = 0*n;                                % initialize calcium
        for j=1:Nc
            C(j:Nc:end) = filter(1,[1, -P.gam(j)],n(j:Nc:end)) + (1-P.gam(j))*P.b(j);
        end

        % precompute parameters required for evaluating and maximizing likelihood
        M(d1)   = -repmat(P.gam,T-1,1);         % matrix transforming calcium into spikes, ie n=M*C
        lam     = dt*repmat(P.lam,T,1);         % for lik
        lnprior = lam.*sum(M)';                 % for grad
        aa      = repmat(diag(P.a'*P.a),T,1);   % for grad
        H1(d0)  = 2*e*aa;                       % for Hess
        gg      = P.a'*F;% for grad
        b       = repmat(P.b,Meta.T,1)';
        bb      = b(:);
        
        % find C = argmin_{C_z} lik + prior + barrier_z
        while z>1e-13                           % this is an arbitrary threshold

            D = F-P.a*(reshape(C,Nc,T)+b);      % difference vector to be used in likelihood computation
            L = e*D(:)'*D(:)+lam'*n-z*sum(log(n));% Likilihood function using C

            s = 1;                              % step size
            d = 1;                              % direction
            while norm(d)>5e-2 && s > 1e-3      % converge for this z (again, these thresholds are arbitrary)
                g   = 2*e*(aa.*(C+bb)-gg(:)) + lnprior - z*M'*(n.^-1);  % gradient
                H2(d0) = n.^-2;                 % part of the Hessian
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
                    D   = F-P.a*(reshape(C1,Nc,T)+b);
                    DD  = D(:)'*D(:);
                    L1  = e*DD+lam'*n-z*sum(log(n));
                    s   = s/10;                  % if step increases objective function, decrease step size
                end
                C = C1;                         % update C
                L = L1;                         % update L
            end
            z=z/10;                             % reduce z (sequence of z reductions is arbitrary)
        end

        % reshape things in the case of multiple neurons within the ROI
        nn=reshape(Z,T,Nc);
        CC=0*nn;
        for k=1:Nc
            nn(:,k) = n(k:Nc:end);
            CC(:,k) = C(k:Nc:end);
        end
        n=nn;C=CC;
    end

end
