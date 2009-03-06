clear, clear global, close all
filehead='Balazs Project#57_T';
filefoot='_B.tif';

nFrames=1000
Cells=zeros(200,nFrames);
for t=2:nFrames
    if t<10
        pad='000';
    elseif t>9 & t<100
        pad='00';
    elseif t>99 & t<1000
        pad='0';
    else
        pad='';
    end
    filename=[filehead pad num2str(t) filefoot];
    folder='C:\D\Classes\CSHL\data\Balazs Project#57.OneFolder\';
    total=[folder filename];
    a=imread(total);
    if t==2
        BWs=edge(a, 'prewitt');
        Ifill=imfill(BWs,'holes');
        se90=strel('line', 3, 90);
        se0=strel('line', 3, 0);
        BWsdil=imdilate(BWs, [se90 se0]);
        BWdfill=imfill(BWsdil,'holes');
        [labeled, numObjects] = bwlabel(BWdfill);

        for n=1:numObjects
            indic{n}=find(labeled==n);
        end

        j=1;
        for n=numObjects:-1:1
            if length(indic{n})>50
                CellNum{j}=indic{n};
                j=j+1;
            end
        end
    end

    for k=1:length(CellNum)
        Cells(k,t)=double(sum(a(CellNum{k})))/length(CellNum{k});
    end
    
    if mod(t,200)==0
        save Cells
        t
    end
end

Cells=Cells(1:197,1:1000);
for n=1:length(CellNum)
    Cells(n,:)=normalize(Cells(n,:),0,1);
end
Burst=round(Cells(1:197,1:1000)*.5);
imagesc(Burst)

%     figure, imshow(BWdfill)
%     BWoutline{t} = bwperim(BWdfill);
%     Segout = a;
%     Segout=uint16(round(normalize(double(Segout),0,2^16-1)));
%     imshow(Segout)
%     Segout(BWoutline{t})=0;
%     figure, imshow(Segout)
% end
%
% Sums=zeros(size(BWdfill{2}));
% for t=2:nFrames
%     Sums=Sums+BWdfill{t};
% end
% Means=Sums/t;
% MeanOutline=round(Means);
% MeanFill=imfill(MeanOutline,'holes');
% [labeled, numObjects] = bwlabel(MeanFill);
% figure, imagesc(labeled); title(numObjects)
%
% [label1, numO1] = bwlabel(BWoutline{2});
% figure, imagesc(label1); title(numO1);
%
% [label2, numO2] = bwlabel(BWdfill{t});
% figure, imagesc(label2); title(numO2);
%
