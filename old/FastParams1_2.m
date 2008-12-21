function P = FastParams1_2(F,C,n,T,dt,tau)
% estimates parameters using mu=0;
%%
P.lam  = sum(n)/(T*dt);

if nargin==5
    W       = C(1:end-1);
    Y       = F(2:end)-n(2:end);
    a       = W\Y;
    if a(1)>1, a(1)=1; elseif a(1)<0, a(1)=0; end           % make sure 'a' is within bounds
    P.tau   = dt/(1-a(1));
elseif nargin==6
    P.tau   = tau;
end

DD      = (F-C)'*(F-C);                            % difference vector
P.sig   = sqrt(DD/T);

c       = 1/(2*P.sig^2);                    % scale of variance
P.lik   = 0.5*T*log(2*pi*P.sig^2) + c*DD - T*log(P.lam*dt) + P.lam*dt*sum(n);% initialize likelihood
P.mu    = 0;

end