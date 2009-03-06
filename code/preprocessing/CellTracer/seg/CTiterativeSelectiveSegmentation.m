function result = CTiterativeSelectiveSegmentation(im,inputmask,para,debug)
cellscore = para{1};
minvolume = para{2};
disk = para{3};
distancepercentile = para{4};
imold = im;
if ~isempty(inputmask)
    im(inputmask(:,:,1) > 0) = 0;
    im(inputmask(:,:,2) > 0) = 0;
    im(inputmask(:,:,3) > 0) = 0;
end
mask = im > 0; 
if disk > 0
   mask = imopen(mask,strel('disk',disk)); 
end
if minvolume > 0
    mask = bwareaopen(mask,minvolume);
    result(mask==0) = 0;
end
mask = imfill(mask,'holes');
im(mask==0) = 0;

L = bwlabel(mask); s = regionprops(L, 'Image','BoundingBox','Solidity','ConvexArea','Area','perimeter');  
for k = 1:numel(s)
    k;
    bb = uint16(s(k).BoundingBox);
    a = im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1); img = a; 
    a(~s(k).Image) = 0; img(a>0) = 0;
    temp = (s(k).ConvexArea - s(k).Area) / s(k).Perimeter;
    flag = temp;
    
    if flag > cellscore
        if debug > 0
            figure(101),imshow(a); title([num2str(flag),'     ', num2str(temp)]); 
            %pause
        end
        b = a;
        abackup =a;
        if minvolume >0
            [dummy newborder]=CTsplit3(b,disk,minvolume);
        elseif minvolume <0
            [dummy newborder]=CTsplit1(b);
        else
            [dummy newborder]=CTsplit2(b);
        end
        b(newborder>0) = 0;
        if disk >0
            smoothed = CTiterativeOpen(b, disk,1);
            smoothed = CTdistanceSmooth(smoothed,disk,distancepercentile,debug);
        else
            smoothed = b;
        end
        b(a>0 & smoothed == 0) = 0;
        %restore
        c = b; c(c>0) = 255; 
        if disk >0
            c = uint8(bwmorph(c,'thicken',disk * 2)); c(c>0) = 255; 
        end
        a(a>0 & c==0) = 0; 
               
        if sum(sum(abackup-a))>0
            a =CTiterativeSelectiveSegmentation(a,[],para,debug);
        end
        
        
        img = img+a;
        im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1) = img;
    end
end

result = im;
if ~isempty(inputmask)
    temp = inputmask(:,:,3);
    if any(temp>0)
        temp = imdilate(temp,strel('disk',1));
        result(temp>0) = 0;
        result(inputmask(:,:,3)>0) = imold(inputmask(:,:,3)>0);
    end
end
if debug > 0
    figure(101),imshow(result)
end
