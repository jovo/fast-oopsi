function [n_best P_best]=FOOPSI_v3_05_01(F,P,Meta,Est)
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
% Input---- only F is REQUIRED.  the others are optional
% F:        fluorescence time series (can be a vector (1 x T) or a matrix (Np x T)
% P.        structure of neuron parameters
%   a:      spatial filter
%   b:      baseline
%   sig:    standard deviation
%   gam:    decayish (ie, tau=dt/(1-gam)
%   lam:    firing rate-ish
% Meta.      structure of simulation parameters
%   T:      # of time steps
%   dt:     time step size
%   Np:     # of pixels in ROI
%   Nc:     # of cells within ROI
%   Poiss:  whether observations are assumed to come from a Poisson or Gaussian distribution
%   MaxIter:maximum number of iterations of pseudo-EM   (typically set to 50)
% Est.      structure of parameters necessary for estimating parameters
%   sig:    whether to estimate sig (default 0)
%   lam:    whether to estimate lam (default 0)
%   gam:    whether to estimate gam (default 0) 
%   b:      whether to estimate b   (default 1)
%   a:      whether to estimate a   (default 1)
%   Plot:   whether to plot results (only required is Est.a==1)
%   Thresh: whether to threshold infered spike train before updating 'a' and 'b' (only required is Est.a==1)
%   n:      if true spike train is known, and we are plotting, plot it (only required is Est.a==1)
%   h:      height of ROI (assumes square ROI) (# of pixels) (only required is Est.a==1)
%   w:      width of ROI (assumes square ROI) (# of pixels) (only required is Est.a==1)
%
% Output---
% n_best:   inferred spike train
% P_best:   inferred parameter structure
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
% 3_03_01:  made background a scalar, inference works, learning does not
% 3_04_01:  added possibility of using Poisson observation noise (but it is
%           still buggy), and made learning work for gaussian observation
% 3_04_02:  modified input structures (see above for details)
% 3_05_01:  lam can be time-varying

%% initialize stuff

% set default parameter values
if nargin < 2, P = struct; end
if ~isfield(P,'b'),     P.b=0; end
if ~isfield(P,'sig'),   P.sig=0.01; end
if ~isfield(P,'gam'),   P.gam=0.975; end
if ~isfield(P,'lam'),   P.lam=10; end
if ~isfield(P,'a'),     [U,S,P.a] = pca_approx(F',1); end

% set meta parameter values
if nargin < 3, Meta = struct; end
if isfield(Meta,'Poiss'),   Poiss  = Meta.Poiss; else Poiss = 0; end
if isfield(Meta,'Nc'),      Nc = Meta.Nc;   else Nc = 1; end
sizF = size(F);
if isfield(Meta,'Np'),      Np = Meta.Np;   else Np = sizF(1); end
if isfield(Meta,'T'),       T  = Meta.T;    else T  = sizF(2); end
if isfield(Meta,'dt'),      dt = Meta.dt;   else
    fr = input('what was the frame rate for this movie (in Hz)? ');
    dt = 1/fr;
end
if isfield(Meta,'MaxIter'), MaxIter = Meta.MaxIter; else
    reply = input('do you want to estimate parameters? y/n [y] (case sensitive): ', 's');
    if reply == 'y'; MaxIter = 10;
    else MaxIter = 0; end
end

% set which parameters to estimate
if nargin < 4, Est = struct; end
if MaxIter>1;
    if ~isfield(Est,'sig'), Est.sig   = 0; end
    if ~isfield(Est,'lam'), Est.lam   = 0; end
    if ~isfield(Est,'gam'), Est.gam   = 0; end
    if ~isfield(Est,'a'),   Est.a     = 1; end
    if ~isfield(Est,'b'),   Est.b     = 1; end
    if isfield(Est,'Thresh'),   Thresh = Est.Thresh; else Thresh = 1; end
    if isfield(Est,'Plot'),     DoPlot = Est.Plot;   else DoPlot = 1; end
else
    Est.a=1;
end

% make sure we have 1 spatial filter per neuron in ROI
siz=size(P.a);
if Est.a==1 && siz(2)~=Nc
    [U,S,V]=pca_approx(F',Nc);
    for j=1:Nc, P.a(:,j)=V(:,j); end
end

%% define some stuff needed for FastFilter function
Z   = zeros(Nc*T,1);                            % zero vector
M   = spdiags([repmat(-P.gam,T,1) repmat(Z,1,Nc-1) (1+Z)], -Nc:0,Nc*T,Nc*T);  % matrix transforming calcium into spikes, ie n=M*C
I   = speye(Nc*T);                              % create out here cuz it must be reused
d0  = 1:Nc*T+1:(Nc*T)^2;                        % index of diagonal elements of TxT matrices
d1  = 1+Nc:Nc*T+1:(Nc*T)*(Nc*(T-1));            % index of off-diagonal elements of TxT matrices
l   = Z(1:MaxIter);                             % initialize likelihood
if numel(P.lam)==Nc
    lam = dt*repmat(P.lam,T,1);         % for lik
elseif numel(P.lam)==Nc*T
    lam = dt*P.lam;
else
    error('lam must either be length T or 1');
end

if Poiss==1
    H       = I;                                % initialize memory for Hessian matrix
    gamlnF  = gammaln(F+1);                     % for lik
    sumF    = sum(F);                           % for Hess
else
    H1  = I;                                    % initialize memory for Hessian matrix
    H2  = I;                                    % initialize memory for Hessian matrix
end

%% infer spike train using default/initialized parameters
n = FastFilter(F,P);

%%  pseudo-EM iterations

% set up stuff
l(1)    = -inf;
l_max   = l(1);                                 % maximum likelihood achieved so far
n_best  = n;                                    % best spike train
P_best  = P;                                    % best parameter estimate
options = optimset('Display','off');
mse_old = -inf;

% iterate until MaxIter or lik stops changing
for i=2:MaxIter

    % update inferred spike train
    [n C]   = FastFilter(F,P);

    if min(n(:))<0,
        disp('somehow, a negative spike has arisen'), keyboard
    end

    % generate regressor for spatial filter
    if Est.a==1 || Est.b==1
        if Thresh==1
            CC=0*C;
            for j=1:Nc
                nsort   = sort(n(:,j));
                nthr    = nsort(round(0.98*T));
                nn      = Z(1:T);
                nn(n(:,j)<=nthr)=0;
                nn(n(:,j)>nthr)=1;
                CC(:,j) = filter(1,[1 -P.gam(j)],nn) + (1-P.gam(j))*P.b(j);
            end
        else
            CC      = C;
        end

        % update spatial filter and baseline
        CC = CC + b';
        if Est.a==1
            for ii=1:Np
                Y   = F(ii,:)';
                P.a(ii,:) = CC\Y;
            end
        end
        if Est.b==1
            P.b     = quadprog(P.a'*P.a,-P.a'*sum(F - P.a*CC',2)/T',[],[],[],[],Z(1:Nc),inf+Z(1:Nc),P.b,options);
            P.b     = P.b';
        end
        b       = repmat(P.b,T,1)';
        D       = F-P.a*(reshape(C,Nc,T)+b);
        mse     = -D(:)'*D(:);
    end

    % estimate other parameters
    if Est.sig==1, P.sig = sqrt(-mse)/T; end
    if Est.lam==1,
        nnorm   = n./repmat(max(n),T,1);
        if numel(P.lam)==Nc
            P.lam   = sum(nnorm)'/(T*dt);
            lam     = repmat(P.lam,T,1)*dt;
        else
            P.lam   = nnorm/(T*dt);
            lam     = P.lam*dt; 
        end
        
    end

    % update likelihood and keep results if they improved
    lik     = -T*Np*log(2*pi*P.sig^2)/2 - mse/(2*P.sig^2);
    prior   = sum(lam(:)) - lam(:)'*n(:);
    l(i)    = lik + prior;

    % if this is the best one, keep n and P
    if l(i)>l_max
        n_best  = n;
        P_best  = P;
        l_max   = l(i);
    end

    % if lik doesn't change much (relatively), or returns to some previous state, stop iterating
    if abs((l(i)-l(i-1))/l(i))<1e-5 || any(l(1:i-1)-l(i))<1e-5% abs((l(i)-l(i-1))/l(i))<1e-5 || l(i-1)-l(i)>1e5;
        disp('modifying theta did not reduce likelihood enough to justify more of it')
        break;
    end

    % plot results from this iteration
    if DoPlot == 1
        if isfield(Est,'h') && isfield(Est,'h')
            figure(400), nrows=Nc;                % plot spatial filter
            for j=1:Nc, subplot(1,nrows,j),
                imagesc(reshape(P.a(:,j),Est.h,Est.w)),
                title('a')
            end
        end

        figure(401), clf, ncols=Nc+1; END=min(T,200);
        for j=1:Nc                              % plot inferred spike train
            h(j)=subplot(ncols,1,j);
            bar(z1(n(2:END,j)))

            if isfield(Est,'n'), hold on,
                stem(Est.n(2:END,j),'LineStyle','none','Marker','v','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',2)
            end
            title(['iteration ' num2str(i)]),
            set(gca,'XTickLabel',[])
            axis('tight')
        end
        subplot(ncols,1,ncols), semilogy(l(2:i))    % plot record of likelihoods
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
        llam = reshape(lam',1,Nc*T)';
        n = z./llam;
        C = 0*n;                                % initialize calcium
        for j=1:Nc
            C(j:Nc:end) = filter(1,[1, -P.gam(j)],n(j:Nc:end)) + (1-P.gam(j))*P.b(j);
        end

        % precompute parameters required for evaluating and maximizing likelihood
        b           = repmat(P.b,T,1)';        % for lik
        if Poiss==1
            suma    = sum(P.a);                     % for grad
        else
            aF      = P.a'*F;                       % for grad
            bb      = b(:);                         % for grad
            M(d1)   = -repmat(P.gam,T-1,1);         % matrix transforming calcium into spikes, ie n=M*C
            lnprior = llam.*sum(M)';                 % for grad
            aa      = repmat(diag(P.a'*P.a),T,1);   % for grad
            H1(d0)  = 2*e*aa;                       % for Hess
        end

        % find C = argmin_{C_z} lik + prior + barrier_z
        while z>1e-13                           % this is an arbitrary threshold

            if Poiss==1
                Fexpected = P.a*(C+b')';        % expected poisson observation rate
                L = sum(sum(exp(-Fexpected+ F.*log(Fexpected) - gamlnF)));
            else
                D = F-P.a*(reshape(C,Nc,T)+b);  % difference vector to be used in likelihood computation
                L = e*D(:)'*D(:)+llam'*n-z*sum(log(n));% Likilihood function using C
            end
            s = 1;                              % step size
            d = 1;                              % direction
            while norm(d)>5e-2 && s > 1e-3      % converge for this z (again, these thresholds are arbitrary)
                if Poiss==1
                    g   = (-suma + sumF./(C+b')')';
                    H(d0) = sumF'.*(C+b').^(-2);
                else
                    g   = 2*e*(aa.*(C+bb)-aF(:)) + lnprior - z*M'*(n.^-1);  % gradient
                    H2(d0) = n.^-2;             % part of the Hessian
                    H   = H1 + z*(M'*H2*M);     % Hessian
                end
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
                    if Poiss==1
                        Fexpected = P.a*(C1+b')';
                        L1 = sum(sum(exp(-Fexpected + F.*log(Fexpected) - gamlnF)));
                    else
                        D   = F-P.a*(reshape(C1,Nc,T)+b);
                        DD  = D(:)'*D(:);
                        L1  = e*DD+llam'*n-z*sum(log(n));
                    end
                    s   = s/5;                  % if step increases objective function, decrease step size
                    if s<1e-20; disp('reducing s further did not increase likelihood'), break; end      % if decreasing step size just doesn't do it
                end
                C = C1;                         % update C
                L = L1;                         % update L
            end
            z=z/10;                             % reduce z (sequence of z reductions is arbitrary)
        end

        % reshape things in the case of multiple neurons within the ROI
        n=reshape(n,Nc,T)';
        C=reshape(C,Nc,T)';
    end

end
