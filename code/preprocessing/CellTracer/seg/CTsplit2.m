function [result,newborder] = CTsplit2(im)
im_backup = im;
D = bwdist(~im);
maxsize = max(D(:));
temp=D(D>0);minsize= min(temp);
split = 1;
while split<2 && minsize < maxsize
    im = im_backup;
    im(D <= minsize) = 0;
    [L split] = bwlabel(im);
    temp = temp(temp>minsize);
    if ~isempty(temp)
        temp = temp(temp>minsize);
        minsize =min(temp);
    else
        minsize = maxsize;
    end
end
if minsize >=maxsize
    result = im_backup;
    newborder = im_backup;newborder(im_backup>0)=0;
else
    % only keep the top 2 sub-blobs
    if split > 2
        a = zeros(split,1);
        for i=1:split
            a(i) = length(find(L==i));
        end
        [sorteda sindex] = sort(a);
        for i=1:split-2
            im(L == sindex(i)) = 0; 
        end
    end

    s1=im;
    temp = im;temp(im>0) = 255;
    temp = uint8(bwmorph(temp,'thicken',uint8(minsize)));
    %figure(102),imshow(temp);

    %deal with the ones that has been left
    s2 = im_backup; s2(temp>0) = 0;

    %figure(103),imshow(s2);

    D = bwdist(~s2);
    s2(D<minsize) = 0;
    temp = s1 + s2;temp(temp>0) = 255;

    temp = uint8(bwmorph(temp,'thicken',uint8(maxsize)));
    temp(im_backup == 0) = 0; 
    result = temp; 
    newborder = im_backup; newborder(result>0) = 0;    
end
