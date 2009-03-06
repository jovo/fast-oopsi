function bordermask = CTborder2(im,mask,para,debug)
maxhalfwidth = para{1};
borderspread = para{2};
disk = para{3};
method = para{4};
if para{5} < 0
    thres = borderspread/2;
else
    thres = para{5};
end
if isempty(mask)
    backgroundmask = zeros(size(im),'uint8');
else
    backgroundmask = uint8(mask(:,:,2) | mask(:,:,1)); 
end
inputmask = ones(size(im),'uint8');
r1 = CTblock(im,inputmask,1-backgroundmask,maxhalfwidth,{@CTranks},{method});
b = uint8(r1 * 255); b(b==0) = 1;
b1 = CTblock(b,inputmask,1-backgroundmask,disk,{@CTcounts},{255 * borderspread,0}); 
if debug > 0
    figure(100),imshow(r1,[]);
    figure(101),imshow(b1,[]);
end
bordermask = b1<thres;
