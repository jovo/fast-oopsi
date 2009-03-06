function bordermask = CTborder3(im,mask,para,debug)
maxhalfwidth = para{1};
borderspread = para{2};
disk = para{3};
minvolume = para{4};
method = para{5};
if para{6} < 0
    thres = borderspread/2;
else
    thres = para{6};
end
if isempty(mask)
    backgroundmask = zeros(size(im),'uint8');
else
    backgroundmask = uint8(mask(:,:,2)>0 | mask(:,:,1) > 0); 
end
outputmask = 1-backgroundmask;
%inputmask = ones(size(im),'uint8');

cum = zeros(size(im),'double'); 
for i = disk:maxhalfwidth
    r1 = CTblock(im,outputmask,outputmask,i,{@CTranks},{method});
    r1(r1<borderspread) = 0;
    temp = r1>0;
    temp = bwareaopen(temp,minvolume);
    r1(temp==0) = 0;
    cum = cum + r1;
    if debug > 0 && max(cum(:)) > 0
        figure(100),imshow(cum/max(cum(:)),[]);
    end
end
r1 = cum/max(cum(:));
% temp = r1 <= thres & backgroundmask == 0;
% temp = imopen(temp,strel('disk',disk));
% bordermask = temp==0 | mask(:,:,1) > 0;
bordermask = r1>thres;
bordermask = bwareaopen(bordermask,minvolume);
if ~isempty(mask)
    bordermask = bordermask > 0 | mask(:,:,1) >0;
end
if debug>0
    figure(101),imshow(bordermask>0);
end
