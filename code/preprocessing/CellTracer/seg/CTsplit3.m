function [result,newborder] = CTsplit3(im,basesize,minvolume)
temp = im>0; 
if sum(double(temp(:))) <= minvolume
    result = im; newborder = im;newborder(im>0)=0;
else
    [dummy InitSplit] = bwlabel(temp,4);
    D = bwdist(~im);
    maxsize = max(D(:));
    temp=D(D>0);minsize= basesize;

    split = InitSplit;
    while split<InitSplit+1 && minsize < maxsize
        im_in = im; im_in(D <= minsize) = 0;
        inc = uint8(minsize);
        if inc == minsize
            inc = inc - 1;
        end
        temp=uint8(bwmorph(im_in,'thicken',inc));
        im_in = im;im_in(temp==0) = 0;
        smoothed = CTiterativeOpen(im_in,basesize,1);
        im_in(smoothed == 0) = 0;
        im_out=im; im_out(im_in>0) = 0;
        
        %make im_out a bit smaller
        smoothed = CTiterativeOpen(im_out,basesize,1);
        im_out(smoothed == 0) = 0;
        bw=im_out;bw(bw>0) = 255;
        perim=bwperim(bw,8);
        im_out(perim>0)=0;

        im_new=im_in + im_out;
        [L n] = bwlabel(im_new,4);

        a = zeros(n,1);
        for i=1:n
            a(i) = length(find(L==i));
            if a(i) < minvolume 
                im_new(L == i) = 0; 
                a(i) = 0;
            end
        end
        split = length(a(a>0)); 
        minsize = minsize+1;
    end
    if minsize >= maxsize || split <=InitSplit
        result = im;
        newborder = im;newborder(im>0)=0;
    else
        temp=uint8(bwmorph(im_new,'thicken',Inf));
        temp = CTiterativeOpen(temp,basesize,1);
        result = im;result(temp==0) = 0;
        newborder = im;newborder(temp>0)=0;

    end
end