% tau=0.1;
% P.gam=1-Meta.dt/tau;
% P.b=0.5;
% n=poissrnd(10*Meta.dt,Meta.T,1);
% C1 = filter(1,[1, -P.gam],n);
% C1 = C1+P.b;
% 
% C2=0*C1; C2(1)=P.b;
% for t=2:Meta.T
%     C2(t)=P.gam*C2(t-1)+Meta.dt*P.b/tau+n(t);
% end
% 
% C3=0*C1; C3(1)=P.b;
% for t=2:Meta.T
%     C3(t)=P.gam*C3(t-1)+(1-P.gam)*P.b+n(t);
% end
% 
% figure(1), clf
% plot(C1,'r')
% hold on
% plot(C2,'--b')
% plot(C3,'-.g')

%%

T=2000;
Nc=2;
Np=300;
F=rand(T,Np);
C1=rand(Nc*T,1);
P.a=rand(Np,Nc);
b=1;
% for t=1:100
%     D   = F-reshape(C1,Nc,T)'*P.a'-b;
% end

for t=1:100
    x   = reshape(C1,Nc,T);
    D   = F-x'*P.a'-b;
end