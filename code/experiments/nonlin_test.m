C=0.01:.1:1e5;
n=1;
kd=180;
% hill=C.^n./(C.^n+kd);
P.n=n; P.k_d=kd;
hill=Hill_v1(P,C);
figure(1), clf, semilogx(hill), hold on

x0(1)=1e-4;
x0(2)=1;
bolt=2./(1+exp(-x0(1)*C+x0(2)));
plot(bolt,'k')

%%

f = @(x)norm(hill-(1./(1+exp(-x(1)*C+x(2)))));
x = fminunc(f,x0);

plot(2./(1+exp(-x(1)*C+x(2)))-1,'g')

%%
x1(1)=.008;
x1(2)=0.1;
plot(2./(1+exp(-x1(1)*C+x1(2))),'r')

%%
% 
% x2(1)=1e-5;
% x2(2)=.1;
% f = @(x)norm(log(hill)-(-x(1)*C+x(2)));
% x = fminunc(f,x2);

x(1)=1e-2;
x(2)=0;
plot(2./(1+exp(-x(1)*C+x(2)))-1,'k')

%%
clf, clear
C=0.01:.1:1e5;
n=1;
kd=180;
a=1;
P.n=n; P.k_d=kd;
hill=Hill_v1(P,C);

semilogx(hill), hold on,
semilogx(-a./(C+kd)+a*C./((C+kd).^2),'k')
semilogx(-a*kd./((C+kd).^2),'g')
