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