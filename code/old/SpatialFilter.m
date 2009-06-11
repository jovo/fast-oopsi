%% this file does everything required to compare a few ways of spatial
%% filters
% 
% 1) extracts a region of interest in a format useful for matlab
% 2) builds a bunch of filters
% 3) plots what the figures look like
% 3) compares the spike extraction accuracy using the different filters

%% 1) extract ROI (this part changes depending on precisely which data set
%% we are getting

clear, clc, %close all
fname='/Users/joshyv/Research/data/rafa/brendon/last_2_weeks/080218a/080218a21.tif';
ROWS=[55 75];
COLS=[65 80];
% MovInf=imfinfo(fname);
for i=1:235 % we get this number from knowing how long the ephys lasted. numel(MovInf)
    X = imread(fname,i,'PixelRegion', {ROWS, COLS});
    DataMat(:,i)=X(:);
end
DataMat=-double(DataMat);


%% 2) make a bunch of filters

MedianFilter=median(DataMat);
MeanFrame=mean(DataMat');
DataMatMeanSub=DataMat-repmat(MeanFrame,size(DataMat,2),1)';
VarFrame=var(DataMat');
[U1,S1,V1] = svd(DataMat,0);
[U,S,V] = svd(DataMatMeanSub,0);

%% 3) plot filters
i=0; 

i=i+1; figure(i), clf,
imagesc(reshape(z1(MeanFrame),diff(ROWS)+1,diff(COLS)+1)),
colormap('gray'); colorbar;

i=i+1; figure(i), clf,
imagesc(reshape(z1(VarFrame),diff(ROWS)+1,diff(COLS)+1)),
colormap('gray'); colorbar;

% i=i+1; figure(i), clf, semilogx(diag(S1))
i=i+1; figure(i), clf,
imagesc(reshape(z1(U1(:,1)),diff(ROWS)+1,diff(COLS)+1)),
colormap('gray'); colorbar;

% i=i+1; figure(i), clf, semilogx(diag(S))
i=i+1; figure(i), clf,
imagesc(reshape(z1(U(:,1)),diff(ROWS)+1,diff(COLS)+1)),
colormap('gray'); colorbar;

%% 4) get fluorescence trace assuming filter

clear R
i=0;

i=i+1;
R{i}.F = MedianFilter;
R{i}.name = 'median';

i=i+1;
R{i}.F = mean(DataMat);
R{i}.name = 'uniform';
R{i}.filter = reshape(z1(MeanFrame),diff(ROWS)+1,diff(COLS)+1);

i=i+1;
R{i}.F = -(DataMat'*MeanFrame')';
R{i}.name = 'mean';

i=i+1;
R{i}.F = (DataMat'*VarFrame')';
R{i}.name = 'var';
R{i}.filter = reshape(z1(VarFrame),diff(ROWS)+1,diff(COLS)+1);

i=i+1;
Denoised = U(:,1)*S(1)*V(:,1)';
R{i}.F = mean(Denoised);
R{i}.name = 'mean subtracted SVD';
R{i}.filter = reshape(z1(U1(:,1)),diff(ROWS)+1,diff(COLS)+1);

i=i+1;
Denoised1 = U1(:,1)*S1(1)*V1(:,1)';
R{i}.F = mean(Denoised1);
R{i}.name = 'SVD';
R{i}.filter = reshape(z1(U(:,1)),diff(ROWS)+1,diff(COLS)+1);

for i=1:length(R)           % normalize all between 0 and 1
    R{i}.F = z1(R{i}.F);
end

%% 5) plot different filtered versions
for j=1:numel(R)
    FF(j,:)=R{j}.F;
end
i=i+1; figure(i), clf, hold on
plot(FF')
legend('median','uniform','mean','var','SVD-mu','SVD')

%% 6) save stuff

load '/Users/joshyv/Research/data/rafa/brendon/last_2_weeks/080218a/D080218a.mat';
j=24;                                               % experiment number

D = D{j};

freq        = 1;
D.spt       = Get_spt(D.V);                         % get spike times
D.n         = SubSampleSpikes(D,freq);              % make spike train
D.data      = DataMat;                              % raw movie
D.T         = size(DataMat,2);                      % # of time steps
D.ROWS      = ROWS;                                 % 
D.COLS      = COLS;                                 

save('FilteredData_D080218a24','R','D')
