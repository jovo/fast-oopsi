function backgroundmask = CTbackground3(im,mask,para,debug)
maxhalfwidth = para{1};
backgroundspread = para{2};
disk = para{3};
greedy = para{4};
fillholes = para{5};

se = strel('disk',maxhalfwidth+1);
erange(:,:,1) = imdilate(im,se);
erange(:,:,2) = imerode(im,se);
bgmask = zeros(size(im),'uint8');
bgmask(im >0 & erange(:,:,1) - erange(:,:,2) < backgroundspread) = 255;
datamask = bgmask == 0;
if fillholes>0
    datamask = imfill(datamask,'holes');
end
datamask = imerode(datamask,strel('disk',maxhalfwidth));
if greedy > 0
    datamask = imerode(datamask,strel('disk',greedy+1));
elseif greedy < 0
    bgmask = datamask ==0;
    bgmask = CTExpand(im,bgmask,4,maxhalfwidth*2,-greedy);
    datamask = bgmask ==0;
end
if disk > 0
    datamask = imopen(datamask,strel('disk',disk));
end
backgroundmask = datamask ==0; %to get the background mask
if ~isempty(mask)
    backgroundmask(mask(:,:,2) > 0) = 1;
    backgroundmask(mask(:,:,1) > 0) = 0;
    backgroundmask(mask(:,:,3) > 0) = 0;
end


