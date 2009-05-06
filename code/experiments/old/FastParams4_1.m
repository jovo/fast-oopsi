function [P U S V] = FastParams4_1(F,Sim)
% this function estimates the parameters of the model:
% F_t = a*C_t + b + sigma*eps_t, eps_t ~ N(0,1)
% C_t = gamma*C_{t-1} + nu + rho*n_t,   n_t ~ Poisson(n_t; p_t)
%
% Inputs
%   F:  fluorescence time series
%   C:  calcium
%   n:  spike train
%   T:  total number of time steps
%   dt: time step size
%   P:  initial parameters, and which ones should we estimate
%
% Outputs
%   P:  updated parameters
%   l:  lik

[U,S,V] = svd(F-repmat(mean(F),Sim.T,1),0);
PCs     = 1:3;
Denoised= U(:,PCs)*S(PCs,PCs)*V(:,PCs)';
Sim.MakMov  = 1;

if Sim.MakMov==1
    G = Denoised;
    G = G-min(G(:));
    G = G/max(G(:));
    G = uint8(round(G*254+1));
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(G(i,:),Sim.h,Sim.w)',[Sim.matname(1:end-4) 'd2.tif'],'tif','Compression','none','WriteMode',mod)
    end
end
figure(5), clf, imagesc(reshape(V(:,1),Sim.h,Sim.w))

P.a     = V(:,1);
P.b     = 0*P.a;

n       = zeros(Sim.T,1); n(Sim.spt)=1;
P.lam   = sum(n)'/(Sim.T*Sim.dt);

% x0      = 0.99;
% Z       = 0*n;
% options = optimset('Display','off');
% for kk=1:10
%     C0      = F*P.a;
%     C0      = C0/max(C0);
%     P.gam   = fmincon(@(x)sum((C0-filter(1,[1 -x],n)).^2),x0,[],[],[],[],0,1,[],options);
C       = filter(1,[1 -Sim.gam],n);
X       = [C 1+Z];
for ii=1:Sim.Np
    Y   = F(:,ii);
    B   = X\Y;
    for j=1:Sim.Nc
        P.a(ii,j) = B(j);
    end
    P.b(ii) = B(end);
end
% end
D       = (F-C*P.a'-(1+Z)*P.b');
P.sig   = sqrt((D(:)'*D(:))/Sim.T);

end