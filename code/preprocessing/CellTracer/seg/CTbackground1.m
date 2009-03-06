function backgroundmask = CTbackground1(im,mask,para,debug)
if para{3} > 0
    [m n] = size(im);
    temp = ones(m+2, n+2,'uint8')*max(im(:));
    temp(2:m+1,2:n+1) = im;
    im = temp;
end
cum = zeros(size(im));
for ballsize = 1:20
    se = strel('disk',ballsize);
    erange(:,:,1) = imdilate(im,se);
    erange(:,:,2) = imerode(im,se);
    temp = double(erange(:,:,1) - erange(:,:,2)) / ballsize;
    cum = cum + temp;
end
cum = (cum - min(cum(:))) / (max(cum(:)) - min(cum(:)));
datamask  = cum>para{1};
if 0 > 1
    datamask = imfill(datamask,'holes');    
else
    datamask = CTiterativeApply(datamask,0,4,{@imfill},'holes');
end
if para{2} > 0
    datamask = imopen(datamask,strel('disk',para{2}));
end
backgroundmask = datamask ==0; %to get the background mask
if ~isempty(mask)
    backgroundmask(mask(:,:,2) > 0) = 1;
    backgroundmask(mask(:,:,1) > 0) = 0;
    backgroundmask(mask(:,:,3) > 0) = 0;
end
%backgroundmask = CTExpand(255-im,backgroundmask,4,20,[]);
if para{3} > 0
    backgroundmask = backgroundmask(2:m+1,2:n+1);
end