function P = FastParams1_3(F,C,n,T,dt,Q)
% estimates parameters using mu=0;
%%
if nargin==6 && isfield(Q,'lam')
    P.lam=Q.lam;
else
    P.lam  = sum(n)/(T*dt);
end

if nargin==6 && isfield(Q,'tau')
    P.tau=Q.tau;
else
    W       = C(1:end-1);
    Y       = F(2:end)-n(2:end);
    a       = W\Y;
    if a(1)>1, a(1)=1; elseif a(1)<0, a(1)=0; end           % make sure 'a' is within bounds
    P.tau   = dt/(1-a(1));
end

DD      = (F-C)'*(F-C);                            % difference vector
if nargin==6 && isfield(Q,'sig')
    P.sig=Q.sig;
else
    P.sig   = sqrt(DD/T);
end

c       = 1/(2*P.sig^2);                    % scale of variance
P.lik   = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood

end