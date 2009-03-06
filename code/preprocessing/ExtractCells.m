clear, clear global, close all

foldname    = 'TS108_6';
dir         = '/Users/joshyv/Research/data/karel/takashi/';
im_dir      = '/Image/';
phys_dir    = '/Physiology/';

for i=3

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
        if j==1, mod='overwrite'; else mod='append'; end
        imwrite(X,[tifname(1:end-4),'b.tif'],'Compression','none','WriteMode',mod)
        DataMat(:,(j+1)/2)=X(:);
    end
end
matname=[tifname(end-14:end-4)];
save(matname,'DataMat','MovInf')
% filehead='Balazs Project#57_T';
% filefoot='_B.tif';
%%
clear
load TS108_6_003
meanframe=reshape(mean(DataMat,2),MovInf(1).Height,MovInf(1).Width);
meanframe=meanframe-min(meanframe(:));
meanframe=floor((2^16-1)*meanframe/max(meanframe(:)));
meanframe=uint16(meanframe);
meantif=[tifname(1:end-4),'_mean.tif'];
imwrite(meanframe,meantif,'Compression','none','WriteMode','overwrite')

%%

[U,S,V] = svd(F,0);
sv1=reshape(V(1,:),MovInf(1).Height,MovInf(1).Width);
sv1=sv1-min(sv1(:));
sv1=floor((2^16-1)*sv1/max(sv1(:)));
sv1=uint16(sv1);
svtif=[tifname(1:end-4),'_sv.tif'];
imwrite(sv1,svtif,'Compression','none','WriteMode','overwrite')

%%
a=imread(svtif);
BWs=edge(a, 'prewitt');
Ifill=imfill(BWs,'holes');
se90=strel('line', 3, 90);
se0=strel('line', 3, 0);
BWsdil=imdilate(BWs, [se90 se0]);
BWdfill{n}=imfill(BWsdil,'holes');

%%
% nFrames=100;
% for n=2:nFrames
%     %     if n<10
%     %         pad='000';
%     %     elseif n>9 & n<100
%     %         pad='00';
%     %     elseif n>99 & n<1000
%     %         pad='0';
%     %     else
%     %         pad='';
%     %     end
%     %     filename=[filehead pad num2str(n) filefoot];
%     %     folder='C:\D\Classes\CSHL\data\Balazs Project#57.OneFolder\';
%     %     total=[folder filename];
%     a=imread(meantif);
%     BWs=edge(a, 'prewitt');
%     %     figure, imshow(BWs)
%     Ifill=imfill(BWs,'holes');
%     %     figure, imshow(Ifill)
%     se90=strel('line', 3, 90);
%     se0=strel('line', 3, 0);
%     BWsdil=imdilate(BWs, [se90 se0]);
%     %     figure; imshow(BWsdil)
%     BWdfill{n}=imfill(BWsdil,'holes');
%     %     figure, imshow(BWdfill)
%     %     BWoutline{n} = bwperim(BWdfill);
%     %     Segout = a;
%     %     Segout=uint16(round(normalize(double(Segout),0,2^16-1)));
%     %     imshow(Segout)
%     %     Segout(BWoutline{n})=0;
%     %     figure, imshow(Segout)
% end
% 
% Sums=zeros(size(BWdfill{2}));
% for n=2:nFrames
%     Sums=Sums+BWdfill{n};
% end
% Means=Sums/n;
% MeanOutline=round(Means);
% MeanFill=imfill(MeanOutline,'holes');
% [labeled, numObjects] = bwlabel(MeanFill);
% figure, imagesc(labeled); title(numObjects)
% 
% [label1, numO1] = bwlabel(BWoutline{2});
% figure, imagesc(label1); title(numO1);
% 
% [label2, numO2] = bwlabel(BWdfill{n});
% figure, imagesc(label2); title(numO2);
% 
