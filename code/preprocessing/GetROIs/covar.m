clc, clear
fact=0.00001;
x=[1:10]+fact*rand(1,10);
y=[1:10]+fact*rand(1,10);

t=-2*pi:.01:2*pi;
x=sin(t);
y=cos(t);
for i=1:length(t)
    covar(i)=(x(i)-mean(x))*(y(i)-mean(y));
end
a=mean(covar)

c=[6:15];
C=repmat(c,5,1)+rand(5,10);
C=C';
siz=size(C)
mn=mean(C,2);
D=C-repmat(mn,1,siz(2));
cova=D*D';
covari=cova/(siz(2)-1);
[PC, V]=eig(covari)

PC(:,end)*V(end)
mnclc 
D


