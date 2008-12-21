function S = GOOPSI_forward_v1_0(Sim,F,P)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that Sim.M=0 (M is the # of spike history terms).
%
% The backward sampler has standard variance, as approximated typically.
% Each particle has the SAME backwards sampler, initialized at E[h_t]
% this function only does spike history stuff if necessary
%
% the model is F_t = f(C) = alpha C^n/(C^n + k_d) + beta + e_t,
% where e_t ~ N[f(C), gamma*f(C)+zeta]
%
% Inputs---
% Sim:  simulation metadata
% F:    fluorescence
% P:    current parameter estimates
%
% Outputs---
% S:    simulation states

%% allocate memory and initialize stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extize particle info
S.p     = zeros(Sim.N,Sim.T);               % extize rate
S.n     = zeros(Sim.N,Sim.T);               % extize spike counts
S.C     = P.C_init*ones(Sim.N,Sim.T);       % extize calcium
S.w_f   = 1/Sim.N*ones(Sim.N,Sim.T);        % extize forward weights
S.w_b   = 1/Sim.N*ones(Sim.N,Sim.T);        % extize forward weights
S.Neff  = 1/Sim.N*ones(1,Sim.T_o);          % extize N_{eff}

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);            % generate intervals
diffs       = ints(2)-ints(1);                  % generate interval size
A.U_resamp  = repmat(ints(1:end-1),Sim.T_o,1)+diffs*rand(Sim.T_o,Sim.N); % resampling matrix

% extize misc stuff
A.U_sampl   = rand(Sim.N,Sim.T);                % random samples
A.epsilon_c = sqrt(P.sig2_c)*randn(Sim.N,Sim.T);% generate noise on C

if Sim.M>0                                      % if spike histories
    S.h = zeros(Sim.N,Sim.T,Sim.M);             % extize spike history terms
    A.epsilon_h = zeros(Sim.N, Sim.T, Sim.M);   % generate noise on h
    for m=1:Sim.M                               % add noise to each h
        A.epsilon_h(:,:,m) = sqrt(P.sig2_h(m))*randn(Sim.N,Sim.T);
    end
else                                            % if not, comput P[n_t] for all t
    S.p = repmat(1-exp(-exp(P.kx)*Sim.dt)',1,Sim.N)';
end

% extize stuff needed for conditional sampling
A.n_sampl   = rand(Sim.N,Sim.T);              % generate random number to use for sampling n_t
A.C_sampl   = rand(Sim.N,Sim.T);              % generate random number to use for sampling C_t
A.oney      = ones(Sim.N,1);                  % vector of ones to call for various purposes to speed things up
A.zeroy     = zeros(Sim.N,1);                 % vector of zeros

% extize stuff needed for REAL backwards sampling
O.p_o       = zeros(2^(Sim.freq-1),Sim.freq);   % extize backwards prob
O.mu_o      = zeros(2^(Sim.freq-1),Sim.freq);   % extize backwards mean
O.sig2_o    = zeros(1,Sim.freq);                % extize backwards variance

% extize stuff needed for APPROX backwards sampling
O.p         = zeros(Sim.freq,Sim.freq);         % extize backwards prob
O.mu        = zeros(Sim.freq,Sim.freq);         % extize backwards mean
O.sig2      = zeros(Sim.freq,Sim.freq);         % extize backwards var

% initialize backwards distributions
s           = Sim.freq;                         % initialize time of next observation
O.p_o(1,s)  = 1;                                % initialize P[F_s | C_s]

[O.mu_o(1,s) O.sig2_o(s)]  = init_lik(P,F(s));

O.p(1,s)    = 1;                                % initialize P[F_s | C_s]
O.mu(1,s)   = O.mu_o(1,s);                      % initialize mean of P[O_s | C_s]
O.sig2(1,s) = O.sig2_o(s);                      % initialize var of P[O_s | C_s]

if Sim.freq>1                                   % if intermitent sampling
    for tt=s:-1:2                               % generate spike binary matrix
        A.spikemat(:,tt-1) = repmat([repmat(0,1,2^(s-tt)) repmat(1,1,2^(s-tt))],1,2^(tt-2))';
    end
    nspikes = sum(A.spikemat')';                % count number of spikes at each time step

    for n=0:Sim.freq-1
        A.ninds{n+1}= find(nspikes==n);         % get index for each number of spikes
        A.lenn(n+1) = length(A.ninds{n+1});     % find how many spikes
    end
end

O = update_moments(Sim,F,P,S,O,A,s);            % recurse back to get P[O_s | C_s] before the first observation

%% do the particle filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=Sim.freq+1:Sim.T-Sim.freq

    if Sim.pf==0 || (F(s+Sim.freq)-P.beta)/P.alpha>0.98 ||(F(s+Sim.freq)-P.beta)/P.alpha<0
        S = prior_sampler(Sim,F,P,S,A,t);       % use prior sampler when gaussian approximation is not good
    else                                        % otherwise use conditional sampler
        S = cond_sampler(Sim,F,P,S,O,A,t,s);
    end

    % at observations
    if mod(t,Sim.freq)==0
        S = strat_resample(Sim,S,t,A.U_resamp); % stratified resample
        O = update_moments(Sim,F,P,S,O,A,t);    % estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
        s = t;                                  % store time of last observation
    end

    if mod(t,100)==0
        if t<1000                               % print # of observations
            fprintf('\b\b\b%d',t)
        elseif t<10000
            fprintf('\b\b\b\b%d',t)
        elseif t<100000
            fprintf('\b\b\b\b\b%d',t)
        end
    end
    
end %for time loop

end %conditional sampler

%% initialize likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu1 sig1] = init_lik(P,F)
%  get the mean (mu1) and variance (sig1) for P[C_t | F_t]
% compute mean
finv    = ((P.k_d*(F-P.beta))./(P.alpha-F+P.beta)).^(1/P.n); %initialize search with f^{-1}(o)
mu1     = finv;
if mu1>0 && imag(mu1)==0
    %     options = optimset('Display','off','GradObj','on','Hessian','on');
    %     mu1     = fminunc(@fnlogL,finv,options);
    sig1=-1/(-(-P.alpha*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)+P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)+2*(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2*(-P.alpha*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)+P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)-(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)*(-P.alpha*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)+P.alpha*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)+3*P.alpha*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2-2*P.alpha*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2-P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)-(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^3*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2+1/2*(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2*(P.gamma*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)-P.gamma*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)-3*P.gamma*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2+2*P.gamma*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2+P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)-1/2*(P.gamma*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)-P.gamma*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)-3*P.gamma*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2+2*P.gamma*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2+P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)+1/2*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2);
else
    mu1=0;
    sig1=0;
end

%     function  [logL dlogL ddlogL] = fnlogL(C)        %this function compute log L = log P(O|H)
%         logL = (((F-fmu_F(C)).^2)./fvar_F(C)+log(fvar_F(C)))/2;
%         if nargout > 1
%             dlogL=-(F-P.alpha*C^P.n/(C^P.n+P.k_d)-P.beta)/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)*(-P.alpha*C^P.n*P.n/C/(C^P.n+P.k_d)+P.alpha*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)+1/2*(F-P.alpha*C^P.n/(C^P.n+P.k_d)-P.beta)^2/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)^2*(2*P.gamma*C^P.n*P.n/C/(C^P.n+P.k_d)-2*P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)-1/2*(P.gamma*C^P.n*P.n/C/(C^P.n+P.k_d)-P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)/(P.gamma*C^P.n/(C^P.n+P.k_d)+P.zeta);
%             if nargout > 2
%                 ddlogL=-(-P.alpha*C^P.n*P.n/C/(C^P.n+P.k_d)+P.alpha*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)^2/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)+2*(F-P.alpha*C^P.n/(C^P.n+P.k_d)-P.beta)/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)^2*(-P.alpha*C^P.n*P.n/C/(C^P.n+P.k_d)+P.alpha*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)*(2*P.gamma*C^P.n*P.n/C/(C^P.n+P.k_d)-2*P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)-(F-P.alpha*C^P.n/(C^P.n+P.k_d)-P.beta)/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)*(-P.alpha*C^P.n*P.n^2/C^2/(C^P.n+P.k_d)+P.alpha*C^P.n*P.n/C^2/(C^P.n+P.k_d)+3*P.alpha*(C^P.n)^2*P.n^2/C^2/(C^P.n+P.k_d)^2-2*P.alpha*(C^P.n)^3/(C^P.n+P.k_d)^3*P.n^2/C^2-P.alpha*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C^2)-(F-P.alpha*C^P.n/(C^P.n+P.k_d)-P.beta)^2/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)^3*(2*P.gamma*C^P.n*P.n/C/(C^P.n+P.k_d)-2*P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)^2+1/2*(F-P.alpha*C^P.n/(C^P.n+P.k_d)-P.beta)^2/(2*P.gamma*C^P.n/(C^P.n+P.k_d)+2*P.zeta)^2*(2*P.gamma*C^P.n*P.n^2/C^2/(C^P.n+P.k_d)-2*P.gamma*C^P.n*P.n/C^2/(C^P.n+P.k_d)-6*P.gamma*(C^P.n)^2*P.n^2/C^2/(C^P.n+P.k_d)^2+4*P.gamma*(C^P.n)^3/(C^P.n+P.k_d)^3*P.n^2/C^2+2*P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C^2)-1/2*(P.gamma*C^P.n*P.n^2/C^2/(C^P.n+P.k_d)-P.gamma*C^P.n*P.n/C^2/(C^P.n+P.k_d)-3*P.gamma*(C^P.n)^2*P.n^2/C^2/(C^P.n+P.k_d)^2+2*P.gamma*(C^P.n)^3/(C^P.n+P.k_d)^3*P.n^2/C^2+P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C^2)/(P.gamma*C^P.n/(C^P.n+P.k_d)+P.zeta)+1/2*(P.gamma*C^P.n*P.n/C/(C^P.n+P.k_d)-P.gamma*(C^P.n)^2/(C^P.n+P.k_d)^2*P.n/C)^2/(P.gamma*C^P.n/(C^P.n+P.k_d)+P.zeta)^2;
%             end
%         end
%     end

%     function mu_F = fmu_F(C)        %this function compute E[F]=f(C)
%         mu_F    = P.alpha*C.^P.n./(C.^P.n+P.k_d)+P.beta;
%     end
% 
%     function var_F = fvar_F(C)      %this function compute V[F]=f(C)
%         var_F   = P.gamma*C.^P.n./(C.^P.n+P.k_d)+P.zeta;
%     end
end %init_lik

%% update moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function O = update_moments(Sim,F,P,S,O,A,t)

s           = Sim.freq;                 % find next observation time

[mu1 sig1]  = init_lik(P,F(t+s));

O.mu_o(1,s) = mu1;                      % initialize mean of P[O_s | C_s]
O.sig2_o(s) = sig1;                     % initialize var of P[O_s | C_s]

O.p(1,s)    = 1;                        % initialize P[F_s | C_s]
O.mu(1,s)   = mu1;                      % initialize mean of P[O_s | C_s]
O.sig2(1,s) = sig1;                     % initialize var of P[O_s | C_s]

if Sim.M>0
    hhat        = zeros(Sim.freq,Sim.M);% extize hhat
    phat        = zeros(1,Sim.freq+1);  % extize phat

    hs          = S.h(:,t,:);           % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)= hs(:,1,1:Sim.M);      % this too
    hhat(1,:)   = sum(repmat(S.w_f(:,t),1,Sim.M).*h,1); % initialize hhat
    phat(1)     = sum(S.w_f(:,t).*S.p(:,t),1);          % initialize phat
end

if Sim.M>0
    for tt=1:s
        % update hhat
        for m=1:Sim.M                   % for each spike history term
            hhat(tt+1,m)=P.g(m)*hhat(tt,m)+phat(tt);
        end
        y_t         = P.kx(tt+t)+P.omega'*hhat(tt+1,:)';% input to neuron
        phat(tt+1)  = 1-exp(-exp(y_t)*Sim.dt);          % update phat
    end
else
    phat  = 1-exp(-exp(P.kx(t+1:t+s)')*Sim.dt); % update phat
end

for tt=s:-1:2
    O.p_o(1:2^(s-tt+1),tt-1)    = repmat(O.p_o(1:2^(s-tt),tt),2,1).*[(1-phat(tt))*ones(1,2^(s-tt)) phat(tt)*ones(1,2^(s-tt))]';
    O.mu_o(1:2^(s-tt+1),tt-1)   = (1-P.a)^(-1)*(repmat(O.mu_o(1:2^(s-tt),tt),2,1)-P.A*A.spikemat(1:2^(s-tt+1),tt-1)-P.a*P.C_0);     %mean of P[O_s | C_k]
    O.sig2_o(tt-1)              = (1-P.a)^(-2)*(P.sig2_c+O.sig2_o(tt)); % var of P[O_s | C_k]

    for n=0:s-tt+1
        nind=A.ninds{n+1};
        O.p(n+1,tt-1)   = sum(O.p_o(nind,tt-1));
        ps          = (O.p_o(nind,tt-1)/O.p(n+1,tt-1))';
        O.mu(n+1,tt-1)  = ps*O.mu_o(nind,tt-1);
        O.sig2(n+1,tt-1)= O.sig2_o(tt-1) + ps*(O.mu_o(nind,tt-1)-repmat(O.mu(n+1,tt-1)',A.lenn(n+1),1)).^2;
    end
end

if s==2
    O.p     = O.p_o;
    O.mu    = O.mu_o;
    O.sig2  = repmat(O.sig2_o,2,1);
    O.sig2(2,2) = 0;
end

while any(isnan(O.mu(:)))              % in case ps=0/0, which yields a NaN, approximate mu and sig
    O.mu(1,:)   = O.mu_o(1,:);
    O.sig2(1,:) = O.sig2_o(1,:);
    ind         = find(isnan(O.mu));
    O.mu(ind)   = O.mu(ind-1)-P.A;
    O.sig2(ind) = O.sig2(ind-1);
end
O.p=O.p+eps; % such that there are no actual zeros
end %function UpdateMoments

%% particle filtering using the prior sampler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = prior_sampler(Sim,F,P,S,A,t)

if Sim.M>0                                      % update noise on h
    for m=1:Sim.M
        S.h(:,t,m)=P.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end

    % update rate and sample spikes
    hs              = S.h(:,t,:);               % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          % this too
    y_t             = P.kx(t)+P.omega'*h';      % input to neuron
    S.p(:,t)        = 1-exp(-exp(y_t)*Sim.dt);  % update rate for those particles with y_t<0
end
S.n(:,t)        = A.U_sampl(:,t)<S.p(:,t);      % sample n
S.C(:,t)        = (1-P.a)*S.C(:,t-1)+P.A*S.n(:,t)+P.a*P.C_0+A.epsilon_c(:,t);% sample C

% get weights at every observation
if mod(t,Sim.freq)==0
    F_mu        = P.alpha*Hill_v1(P,S.C(:,t))+P.beta;   % compute E[F_t]
    F_var       = P.gamma*Hill_v1(P,S.C(:,t))+P.zeta;   % compute V[F_t]
    ln_w        = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2;% compute log of weights
    ln_w        = ln_w-max(ln_w);                       % subtract the max to avoid rounding errors
    w           = exp(ln_w);                            % exponentiate to get actual weights
    S.w_f(:,t)  = w/sum(w);                             % normalize to define a legitimate distribution
end

end

%% particle filtering using the CONDITIONAL sampler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = cond_sampler(Sim,F,P,S,O,A,t,s)

% if spike histories, sample h and update p
if Sim.M>0                                  % update noise on h
    for m=1:Sim.M                           % for each spike history term
        S.h(:,t,m)  = P.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end
    hs              = S.h(:,t,:);           % this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);      % this too
    S.p(:,t)        = 1-exp(-exp(P.kx(t)+P.omega'*h')*Sim.dt);% update p
end

% compute P[n_k | h_k]
ln_n    = [log(S.p(:,t)) log(1-S.p(:,t))];  % compute [log(spike) log(no spike)]

% compute log G_n(n_k | O_s) for n_k=1 and n_k=0
k   = Sim.freq-(t-s)+1;

m0  = (1-P.a)*S.C(:,t-1)+P.a*P.C_0;         % mean of P[C_k | C_{t-1}, n_k=0]
m1  = (1-P.a)*S.C(:,t-1)+P.A+P.a*P.C_0;     % mean of P[C_k | C_{t-1}, n_k=1]

m2  = O.mu(1:k,t-s);                        % mean of P[O_s | C_k] for n_k=1 and n_k=0
v2  = O.sig2(1:k,t-s);                      % var of P[O_s | C_k] for n_k=1 and n_k=0
v   = repmat(P.sig2_c+v2',Sim.N,1);         % var of G_n(n_k | O_s) for n_k=1 and n_k=0

ln_G0= -0.5*log(2*pi.*v)-.5*(repmat(m0,1,k)-repmat(m2',Sim.N,1)).^2./v;   % log G_n(n_k | O_s) for n_k=1 and n_k=0
ln_G1= -0.5*log(2*pi.*v)-.5*(repmat(m1,1,k)-repmat(m2',Sim.N,1)).^2./v;   % log G_n(n_k | O_s) for n_k=1 and n_k=0

mx  = max(max(ln_G0,[],2),max(ln_G1,[],2))';% get max of these
mx  = repmat(mx,k,1)';

G1  = exp(ln_G1-mx);                        % norm dist'n for n=1;
M1  = G1*O.p(1:k,t-s);                      % times prob of n=1

G0  = exp(ln_G0-mx);                        % norm dist'n for n=0;
M0  = G0*O.p(1:k,t-s);                      % times prob n=0

ln_G    = [log(M1) log(M0)];                % ok, now we actually have the gaussians

% compute q(n_k | h_k, O_s)
ln_q_n  = ln_n + ln_G;                      % log of sampling dist
mx      = max(ln_q_n,[],2);                 % find max of each column
mx2     = repmat(mx,1,2);                   % matricize
q_n     = exp(ln_q_n-mx2);                  % subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
q_n     = q_n./repmat(sum(q_n,2),1,2);      % normalize to make this a true sampling distribution (ie, sums to 1)

% sample n
S.n(:,t)= A.n_sampl(:,t)<q_n(:,1);          % sample n
sp      = S.n(:,t)==1;                      % store index of which samples spiked
nosp    = S.n(:,t)==0;                      % and which did not

% sample C
if mod(t,Sim.freq)==0                       % if not intermittent
    v       = repmat(O.sig2(1,t-s),Sim.N,1);% get var
    m       = repmat(O.mu(1,t-s),Sim.N,1);  % get mean
else                                        % if intermittent, sample from mixture
                                            % first sample component
    [fo,sp_i]   = histc(A.C_sampl(sp,t),[0  cumsum(O.p(1:k-1,t-s))'/sum(O.p(1:k-1,t-s))]);
    [fo,nosp_i] = histc(A.C_sampl(nosp,t),[0  cumsum(O.p(1:k,t-s))'/sum(O.p(1:k,t-s))]);

    v       = O.sig2(1:k,t-s);              % get var of each component
    v(sp)   = v(sp_i);                      % if particle spiked, then use variance of spiking
    v(nosp) = v(nosp_i);                    % o.w., use non-spiking variance

    m       = O.mu(1:k,t-s);                % get mean of each component
    m(sp)   = m(sp_i);                      % if particle spiked, use spiking mean
    m(nosp) = m(nosp_i);                    % o.w., use non-spiking mean
end
v_c         = (1./v+1/P.sig2_c).^(-1);      %variance of dist'n for sampling C
m_c         = v_c.*(m./v+((1-P.a)*S.C(:,t-1)+P.A*S.n(:,t)+P.a*P.C_0)/P.sig2_c);%mean of dist'n for sampling C
S.C(:,t)    = normrnd(m_c,sqrt(v_c));       % sample C

% update weights
if mod(t,Sim.freq)==0                       % at observations compute P(O|H)
    F_mu        = P.alpha*Hill_v1(P,S.C(:,t))+P.beta;           % compute E[F_t]
    F_var       = P.gamma*Hill_v1(P,S.C(:,t))+P.zeta;           % compute V[F_t]
    log_PO_H    = -0.5*(F(t)-F_mu).^2./F_var - log(F_var)/2;    % compute log of weights
else                                        % when no observations are present, P[O|H^{(i)}] are all equal
    log_PO_H    = (1/Sim.N)*A.oney;
end
log_n           = A.oney;                   % extize log sampling spikes
log_n(sp)       = log(S.p(sp,t));           % compute log P(spike)
log_n(nosp)     = log(1-S.p(nosp,t));       % compute log P(no spike)
log_C_Cn        = -0.5*(S.C(:,t)-((1-P.a)*S.C(:,t-1)+P.A*S.n(:,t)+P.a*P.C_0)).^2/P.sig2_c;%log P[C_k | C_{t-1}, n_k]

log_q_n         = A.oney;                   % initialize log q_n
log_q_n(sp)     = log(q_n(sp,1));           % compute what was the log prob of sampling a spike
log_q_n(nosp)   = log(1-q_n(nosp,1));       % or sampling no spike
log_q_C         = -0.5*(S.C(:,t)-m_c).^2./v_c;  % log prov of sampling the C_k that was sampled

log_quotient    = log_PO_H + log_n + log_C_Cn - log_q_n - log_q_C;

sum_logs        = log_quotient+log(S.w_f(:,t-1));   % update log(weights)
w               = exp(sum_logs-max(sum_logs));      % exponentiate log(weights)
S.w_f(:,t)      = w./sum(w);                        % normalize such that they sum to unity

end %condtional sampler

%% stratified resample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = strat_resample(Sim,S,t,U_resamp)
% Inputs:
% Sim       : simulation parameters (eg, dt, N, etc.)
% S         : simulation states (eg, n, lambda, C, h)
% t         : current time step index
% U_resamp  : resampling matrix (so that i needed generate random numbers
%             which each resample, but rather, can just call them
% Outputs a structure S, which has particle states having resampled them, and equalized weights

Nresamp=t/Sim.freq;                         % increase sample counter
S.Neff(Nresamp)  = 1/sum(S.w_f(:,t).^2);    % store N_{eff}

% if weights are degenerate or we are doing prior sampling then resample
if S.Neff(Nresamp) < Sim.N/2 || Sim.pf==0
    [foo,ind]   = histc(U_resamp(Nresamp,:),[0  cumsum(S.w_f(:,t))']);
    [ri,ri]     = sort(rand(Sim.N,1));      % these 3 lines stratified resample
    ind         = ind(ri);

    S.p(:,t-Sim.freq+1:t)   = S.p(ind,t-Sim.freq+1:t);      % resample probabilities (necessary?)
    S.n(:,t-Sim.freq+1:t)   = S.n(ind,t-Sim.freq+1:t);      % resample calcium
    S.C(:,t-Sim.freq+1:t)   = S.C(ind,t-Sim.freq+1:t);      % resample calcium
    S.w_f(:,t-Sim.freq+1:t) = 1/Sim.N*ones(Sim.N,Sim.freq); % reset weights
    if Sim.M>0                                              % if spike history terms
        S.h(:,t-Sim.freq+1:t,:) = S.h(ind,t-Sim.freq+1:t,:);% resample all h's
    end
end %function
end