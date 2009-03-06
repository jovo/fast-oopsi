function [result,newborder] = CTsplit1(im)
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
    if length(temp) > 0
        temp = temp(temp>minsize);
        minsize =min(temp);
    else
        minsize = maxsize;
    end
end
im = uint8(bwmorph(im,'thicken',Inf));
im(im_backup == 0) = 0; 
result = im; 
newborder = im_backup; newborder(im>0) = 0;