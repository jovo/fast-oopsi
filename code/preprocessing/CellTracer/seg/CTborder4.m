function result = CTborder4(im,inputmask,para,debug)
maxsize = para{1};
minpercentile = para{2};
disk = para{3};
maxpercentile = para{4};
globalthreshold = para{5};

imold = im;
if ~isempty(inputmask)
    im(inputmask(:,:,1) > 0) = 0;
    im(inputmask(:,:,2) > 0) = 0;
    im(inputmask(:,:,3) > 0) = 0;
end

%absolute global thresholding
if globalthreshold <255
    im(im>=globalthreshold) = 0;
end

%do local thresholding only if a local threshold is defined
if maxsize > 0
    mask = im > 0; L = bwlabel(mask,4); s = regionprops(L, 'Image','BoundingBox');  
    for k = 1:numel(s)
        bb = uint16(s(k).BoundingBox);
        a = im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1); img = a; 
        
        a(~s(k).Image) = 0; img(a>0) = 0;
        
        b = a;
        localmask = a>0;
        estat = CTblock(b,localmask,localmask,maxsize, {@mean},[]);
        Mean = double(b) - estat(:,:,1);
        Mean = Mean - min(Mean(:));
        r1 = CTblock(Mean,localmask,localmask,maxsize,{@CTranks},{1});
        localmask(r1 > maxpercentile) = 0;
        a(localmask==0) = 0;
        b(localmask==0) = 0;
        temp = r1 < minpercentile & localmask > 0;
        if disk >0
            smoothed = CTiterativeOpen(temp, disk,1);
%             smoothed = CTdistanceSmooth(smoothed,0,0,debug);
%             if disk > 2
%                 smoothed = CTiterativeOpen(smoothed, 1,0);
%             end
        else
            smoothed = temp;
        end
        b(a>0 & smoothed == 0) = 0;
        
        %restore
        c = b>0;  c = uint8(bwmorph(c,'thicken',maxsize)); 
        a(a>0 & c==0) = 0;  
        
        if disk >0
            temp = a>0;
            smoothed = CTiterativeOpen(temp, disk,1);
            a(smoothed==0) = 0;
        end
        img = img+a; 
        im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1) = img;
    end
end
result = imold;
result(im>0) = 0; result(result>0) = 255;
if ~isempty(inputmask)
    result(inputmask(:,:,1) > 0) = 255;
end
