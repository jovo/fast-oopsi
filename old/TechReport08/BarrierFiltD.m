function [n P]=BarrierFiltD(F,E,Sim)
% this function solves the following optimization problem:
% n = argmin sum_t ((F_t - C_t)^2 + lambda n_t
% with constraints: n_t>0 forall t, and
% C_t = a C_{t-1} + A n_t, where a=(1-dt/tau), and A is fixed at 1
%
% the approach is called the barrier method, or log barrier method.
% essentially, we iteratively solve
%
% n_eta = argmin sum_t ((F_t - C_t)^2 + lambda n_t -eta*log(n_t)
%
% starting with eta at some reasonably large value, and then decreasing.
%
% Input----
% F:    fluorescence time series
% E.    structure of neuron parameters
%   tau: time constant size
%   lam: prior weight = 1/(rate*A*dt)
%   sig: standard deviation of err
% Sim.  structure of simulation parameters
%   dt      : time step size
%   T       : number of time steps
%   Mstep   : whether to estimate parameters
%   Plot    : whether to plot results
%
% Output---
% n:    inferred spike train
% P:    inferred parameter structure
%
% note that in this version we assume that we know lambda, so we do not
% estimate it.

dt  = Sim.dt;
T   = Sim.T;
a   = 1 - dt/E.tau;                         %decay factor
o   = 1+0*F;                                %init a unity vector
M   = spdiags([-a*o o], -1:0,T,T);          %matrix transforming calcium into spikes, ie n=M*C
C   = o;                                    %initialize calcium guess
N   = M*C;                                  %initialize spike train
Tmat=speye(T);                              %create out here cuz it must be reused
c   = 1/(2*E.sig^2);                        %scale of variance
Hmat= 2*c*Tmat;                             %Hessian

lik     = zeros(1,50);                      %initialize Likilihood function
lik(1)  = c*(F-C)'*(F-C)+E.lam*sum(N);
i       = 2;
minlik  = inf;
conv    = false;                            %param est converge?
while conv == false

    eta = 1;                                %weight on barrier function
    C   = o;                                %initialize calcium guess

    while eta>1e-13
        D = F-C;                            %difference vector
        N = M*C;                            %get spike train assuming a particular C
        L = c*D'*D+E.lam*sum(N)-eta*sum(log(N));  %Likilihood function using C
        oldL = Inf;
        s = 1;                              %step size
        d = 1;                              %direction
        while s > 1e-3 && norm(d)>5e-2 && L-oldL<0%converge for this eta
            g   = -2*c*D+E.lam*sum(M)'-eta*M'*(N.^-1);%gradient
            H   = Hmat+2*eta*M'*spdiags(N.^-2,0,T,T)*M;%Hessian
            d   = -H\g;                     %direction to step using newton-raphson
            hit = -(N)./(M*d);              %step within constraint boundaries
            hit(hit<0)=[];                  %ignore negative hits
            if any(hit<1)
                s = min(1,0.99*min(hit(hit>0)));
            else
                s = 1;
            end
            L_new = L+1;
            while L_new-L>=1e-2             %make sure newton step doesn't increase objective
                C_new   = C+s*d;
                N       = M*C_new;
                D       = F-C_new;
                L_new   = c*D'*D+E.lam*sum(N)-eta*sum(log(N));
                s       = s/2;              %if step increases objective function, decrease step size
            end
            if any(M*(C+2*s*d)<0), keyboard, end
            C = C+2*s*d;
            N = M*C;

            D = F-C;
            oldL = L;
            L = c*D'*D+E.lam*sum(N)-eta*sum(log(N));
        end
        eta=eta/10;                         %reduce eta
    end

    if Sim.Mstep==true;
        lik(i)     = c*(F-C)'*(F-C)+E.lam*sum(N);
        if i<50 %lik(i)-lik(i-1)<=-1e-4
%             if Sim.Plot == true
%                 figure(99)
%                 subplot(211), hold on, plot(i,lik(i),'o'), axis('tight')
%                 subplot(212), cla, hold on, plot(F,'k'), plot(C,'b'), drawnow
%             end
            
            if lik(i)<minlik, minlik=lik(i); P = E; n = N; end

            i = i+1;
            W = -C(1:end-1);
            Y = -F(2:end)+n(2:end);
            a = W'*Y/(W'*W);
            if a>1, a=1; elseif a<0, a=0; end   %make sure 'a' is within bounds

            E.tau = Sim.dt/(1-a);
            E.sig = sqrt((F-C)'*(F-C)/T);
            %E.lam = T*Sim.dt/sum(n);

            c     = 1/(2*E.sig^2);
            Hmat  = 2*c*Tmat;
            M     = spdiags([-a*o o], -1:0,T,T);          %matrix transforming calcium into spikes, ie n=M*C
        else
            conv = true;
        end
    else
        n = N;
        P = E;
        conv = true;
    end
end