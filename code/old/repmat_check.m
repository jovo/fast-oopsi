mu1 = 1;
mu2 = 2;
sig1 = 1;
sig2 = 2;
OO = ones(1000,1);
x1 = randn(mu1*OO,sig1*OO);
x2 = randn(mu2*OO,sig2*OO);

a = .8;

x=a*x1+(1-a)*x2;

hist(x)


%%

ass=rand(10000,1);
OO=ones(1,100);
clc
tic
for j=1:100
    b=repmat(ass,1,100);
end
toc
tic
for j=1:100
    b=ass*OO;
end
toc