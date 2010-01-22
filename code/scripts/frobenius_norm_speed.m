clear, clc

T=5000;
N=5000;
J=5;
A=rand(N,1);
C=rand(1,T);
e=0.1*rand(N,T);
F=A*C+e;
time1=zeros(J,1);
time2=zeros(J,1);

a=repmat(A.^2,1,T);
g=F./repmat(A,1,T);
CC=repmat(C,N,1);
for j=1:J
    starttime=cputime;
    mse1=sum(sum((g-CC).^2.*a));
    time2(j)=cputime-starttime;

    starttime=cputime;
    d0=F-A*C;
    mse0=d0(:)'*d0(:);
    time1(j)=cputime-starttime;

    err(j)=norm(mse0-mse1);
end

time1, time2
% mse0=0;
% for t=1:T
%     for n=1:N
%         mse0 = mse0 + (F(n,t)-A(n)*C(t))^2;
%     end
% end
% mse0
