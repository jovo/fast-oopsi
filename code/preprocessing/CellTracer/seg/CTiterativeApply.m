function result = CTiterativeApply(im,ballsize,connection,method,para)
mask = im > 0; L = bwlabel(mask,connection); s = regionprops(L, 'Image','BoundingBox');  
for k = 1:numel(s)
    bb = uint32(s(k).BoundingBox);
    a = im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1); img = a; 
    a(~s(k).Image) = 0; img(a>0) = 0;
    %figure(101),imshow(a);
    
    b = a;
    if ballsize > 0
        t = CTblock(b,[],[],ballsize,method,para);
        t = uint8(t * 255);
        t(t==0 & b>0) = 1;
    else
        if ~isempty(para)
            t = method{1}(b,para);
        else
            t = method{1}(b);
        end
    end
    img = img+t;
    im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1) = img;
end

result = im;