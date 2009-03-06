clear, clc, %close all
foldname    = 'TS108_6';
dir         = '/Users/joshyv/Research/data/karel/takashi/';
im_dir      = '/Image/';
phys_dir    = '/Physiology/';

i=9

% get whole movie
if i<10                                         % get tif file name
    tifname=[dir foldname im_dir  foldname '_00' num2str(i) '.tif'];
else
    tifname=[dir foldname im_dir  foldname '_0' num2str(i) '.tif'];
end
MovInf  = imfinfo(tifname);                     % get number of frames
T       = numel(MovInf)/2;                      % only alternate frames have functional data
DataMat = zeros(MovInf(1).Width*MovInf(1).Height,T);% initialize mat to store movie
for j=1:2:T*2
    X = imread(tifname,j);%,'PixelRegion'), {ROWS, COLS});
    imwrite(X,'ass','tif','Compression','none','WriteMode','append')
%     DataMat(:,(j+1)/2)=X(:);
end


%%
PCs=1:3; 
Denoised = U(:,PCs)*S(PCs,PCs)*V(:,PCs)';
plot(mean(Denoised))
R{i}.F = mean(Denoised);
figure(1000), clf, plot(R{i}.F)

XX=Denoised;
XX=XX-min(XX(:));
XX=floor(255*XX/max(XX(:)));
figure(1001), hist(XX(:),100)
XX=uint8(XX);
figure(1002), clf, imagesc(XX)
i=i+1;
for j=1:500
    imwrite(reshape(XX(:,j),Nrows,Ncols)',['ass' num2str(i) '.tif'],'tif','Compression','none','WriteMode','append')
end


%%