function Z = GOOPSI_backward_v1_0(Sim,S,P,Z,t)
% this function iterates backward one step computing P[H_t | H_{t+1},O_{0:T}] 
% Input---
% Sim:  simulation metadata
% S:    particle positions and weights
% P:    parameters
% Z:    a bunch of stuff initialized for speed
% t:    current time step
% 
% Output is a single structure Z with the following fields
% n1:   vector of spikes or no spike for each particle at time t
% C0:   calcium positions at t-1
% C1:   calcium positions at t (technically, this need not be output)
% C1mat:matrix from C1
% C0mat:matrix from C0
% w_b:  backwards weights

% compute ln P[n_t^i | h_t^i]
Z.n1            = S.n(:,t);                         %for prettiness sake
ln_Pn           = 0*Z.oney;                         %for fastiness sake
ln_Pn(Z.n1==1)  = log(S.p(Z.n1==1,t));              %P[n=1] for those that spiked
ln_Pn(~Z.n1)    = log(1-S.p(~Z.n1,t));              %P[n=0] for those that did not

% compute ln P[C_t^i | C_{t-1}^j, n_t^i]
Z.C0        = S.C(:,t-1);                           %for prettiness sake
Z.C1        = S.C(:,t);
Z.C1mat     = Z.C1(:,Z.oney);                       %recall from previous time step
Z.C0mat     = Z.C0(:,Z.oney);                       %faster than repamt
mu          = (1-P.a)*S.C(:,t-1)+P.A*Z.n1+P.a*P.C_0;%mean
mumat       = mu(:,Z.oney)';                        %faster than repmat
ln_PC_Cn    = -0.5*(Z.C1mat - mumat).^2/P.sig2_c;   %P[C_t^i | C_{t-1}^j, n_t^i]

% compute ln P[h_t^i | h_{t-1}^j, n_{t-1}^i]
ln_Ph_hn    = Z.zeroy;                              %reset transition prob for h terms
for m=1:Sim.M                                       %for each h term
    h1      = S.h(:,t,m);                           %faster than repmat
    h1      = h1(:,Z.oney);
    h0      = P.g(m)*S.h(:,t-1,m)+S.n(:,t-1);
    h0      = h0(:,Z.oney)';
    ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/P.sig2_h(m);
end

% compute P[H_t^i | H_{t-1}^j]
sum_lns = ln_Pn(:,Z.oney)+ln_PC_Cn + ln_Ph_hn;      %in order to ensure this product doesn't have numerical errors
mx      = max(sum_lns,[],1);                        %find max in each of row
mx      = mx(Z.oney,:);                             %make a matrix of maxes
T0      = exp(sum_lns-mx);                          %exponentiate subtracting maxes (so that in each row, the max entry is exp(0)=1
Tn      = sum(T0,1);                                %then normalize
T       = T0.*repmat(1./Tn(:)', Sim.N, 1);          %such that each column sums to 1

% compute P[H_t^i, H_{t-1}^j | O]
PHHn    = (T*S.w_f(:,t-1))';                        %denominator
PHHn(PHHn==0) = eps;
PHHn2   = PHHn(Z.oney,:)';                          %faster than repmat
PHH     = T .* (S.w_b(:,t)*S.w_f(:,t-1)')./PHHn2;   %normalize such that sum(PHH)=1
sumPHH  = sum(PHH(:));
if sumPHH==0
    Z.PHH = ones(Sim.N)/(Sim.N);
else
    Z.PHH   =  PHH/sum(PHH(:));
end
Z.w_b   = sum(Z.PHH,1);                             %marginalize to get P[H_t^i | O]

if any(isnan(Z.w_b))
    return
end

if mod(t,100)==0 && t>=9900
    fprintf('\b\b\b\b\b%d',t)
elseif mod(t,100)==0 && t>=900
    fprintf('\b\b\b\b%d',t)
elseif mod(t,100)==0
    fprintf('\b\b\b%d',t)
end
