function post = GetLik2_0(DD,n,P,Sim)
% se=zeros(Sim.Np,Sim.T);
% for ii=1:Sim.Np
%     for t=1:Sim.T
%         se(ii,t)=(F(t,ii)-P.a(ii,:)*C(t,:)'-P.b(ii))^2;
%     end
% end
sound(3*sin(linspace(0,90*pi,2000)))
% lik     = -Sim.T*Sim.Np*log(2*pi*P.sig^2)/2 -1/(2*P.sig^2)*sum(se(:));
lik     = -Sim.T*Sim.Np*log(2*pi*P.sig^2)/2 -1/(2*P.sig^2)*DD;
prior   = Sim.T*sum(P.lam*Sim.dt) - Sim.dt*P.lam'*sum(n)';
post    = lik + prior;

%     aC=P.a*C';
%     b=repmat(P.b,1,Sim.T);
%     mse2=F(:)-aC(:)-b(:);
%     sum(mse2.^2)
%     mse3=D(:)'*D(:)
%     0.5*T*log(2*pi*P.sig^2) + norm(F'-P.a*C'-P.b*Z(1:T)')^2/(2*P.sig^2)- T*log(P.lam*dt) + P.lam*dt*sum(n);
