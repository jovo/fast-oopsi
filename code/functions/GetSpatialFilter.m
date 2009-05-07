function R = GetSpatialFilter(F,Nrows,Ncols)

%% 1) make a bunch of filters

MedianFilter=median(F);
MeanFrame=mean(F');
DataMatMeanSub=F-repmat(MeanFrame,size(F,2),1)';
VarFrame=var(F');
[U1,S1,V1] = svd(F,0);
[U,S,V] = svd(DataMatMeanSub,0);

%% 2) plot filters
i=0; 

i=i+1; figure(i), clf,
imagesc(reshape(z1(MeanFrame),Nrows,Ncols)'),
colormap('gray'); colorbar;

i=i+1; figure(i), clf,
imagesc(reshape(z1(VarFrame),Nrows,Ncols)'),
colormap('gray'); colorbar;

% i=i+1; figure(i), clf, semilogx(diag(S1))
i=i+1; figure(i), clf,
imagesc(reshape(z1(U1(:,1)),Nrows,Ncols)'),
colormap('gray'); colorbar;

% i=i+1; figure(i), clf, semilogx(diag(S))
i=i+1; figure(i), clf,
imagesc(reshape(z1(U(:,1)),Nrows,Ncols)'),
colormap('gray'); colorbar;

%% 3) get fluorescence trace assuming filter

clear R
i=0;

i=i+1;
R{i}.F = MedianFilter;
R{i}.name = 'median';

i=i+1;
R{i}.F = mean(F);
R{i}.name = 'uniform';
R{i}.filter = reshape(z1(MeanFrame),Nrows,Ncols);

i=i+1;
R{i}.F = -(F'*MeanFrame')';
R{i}.name = 'mean';

i=i+1;
R{i}.F = (F'*VarFrame')';
R{i}.name = 'var';
R{i}.filter = reshape(z1(VarFrame),Nrows,Ncols);

i=i+1;
Denoised = U(:,1)*S(1)*V(:,1)';
R{i}.F = mean(Denoised);
R{i}.name = 'SVD';
R{i}.filter = reshape(z1(U1(:,1)),Nrows,Ncols);

i=i+1;
Denoised1 = U1(:,1)*S1(1)*V1(:,1)';
R{i}.F = mean(Denoised1);
R{i}.name = 'mean subtracted SVD';
R{i}.filter = reshape(z1(U(:,1)),Nrows,Ncols);

for i=1:length(R)           % normalize all between 0 and 1
    R{i}.F = -z1(R{i}.F)+1;
end

%% 4) plot different filtered versions
for j=1:numel(R)
    FF(j,:)=R{j}.F;
end
i=i+1; figure(i), clf, hold on
plot(FF')
% legend('median','uniform','mean','var','SVD-mu','SVD')