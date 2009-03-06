function rgb = drawimage(currentview, whichaxis,im,ccolor,labels)
    if currentview<3 || isempty(ccolor) % not tracking view
        if max(im(:)) ==1
            rgb = uint8(im>0) * 255;
        else
            rgb = im;
        end
    else
        im = uint8(im) * 255;
        t1 = im; t2 = im; t3 = im;
        curerntcolors = ccolor{whichaxis};
        csize = size(curerntcolors);
        for i = 1:csize(1);
            t1(labels==curerntcolors(i,1)) = curerntcolors(i,2);
            t2(labels==curerntcolors(i,1)) = curerntcolors(i,3);
            t3(labels==curerntcolors(i,1)) = curerntcolors(i,4);
        end
        rgb(:,:,1) = t1; rgb(:,:,2) = t2; rgb(:,:,3) = t3;
    end
end