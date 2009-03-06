clc
clear
close all

vector=[1:.5:60];
T=length(vector);
N=10;
data=repmat(vector, N, 1) + 2*rand(N,T);
plot(data,vector)
mn=mean(data,2);
data=data-repmat(mn,1,T);

covariance1 = 1/(N-1) * data * data';

[PC1 V1]=eig(covariance1);

v1=diag(V1);
[junk, order1] = sort(-1*v1);

v1=v1(order1);
PC1 = PC1(:,order1);
V1=diag(v1,0);

xinv=(sqrt(V1)*PC1')^-1;
PC3=data'*xinv;

reconstruct1=PC1*V1*PC3';
reconstruct1=normalize(reconstruct1,1,60);
figure, plot(reconstruct1,vector)

covariance2 = 1/(N-1) * data' * data;
[PC2 V2]=eig(covariance2);

v2=real(diag(V2));
[junk, order2] = sort(-1*v2);

v2=v2(order2);
PC2 = PC2(:,order2);
V2=diag(v2,0);

xinv2=(sqrt(V2)*PC2')^-1;
PC4=data*xinv2;

reconstruct2=real(PC2*V2*PC4');
reconstruct2=normalize(reconstruct2,1,60)';
figure, plot(reconstruct2,vector)

data=normalize(data,1,60);
figure, plot(data,vector)


figure,
for i=1:size(PC3,2)
    subplot(2,5,i), plot(vector,PC3(:,i)), axis([0 60 -.8 .8])
end
    