function [t1 t2] = CTshowGlobalAlignment(im1,im2,alignment,fig)
[m n] = size(alignment.score); shift = (m+1)/2;
bestlocation = alignment.bestlocation;
[m n] = size(im1);container = zeros(m+2* shift, n+ 2 * shift);
t1 = container; t1(shift+1:shift+m, shift+1:shift+n) = im1;
t2 = container; t2(bestlocation(1):bestlocation(1) + m-1, bestlocation(2):bestlocation(2)+n-1) = im2;
if fig>0
    rgb(:,:,1) = t1;
    rgb(:,:,2) = t2;
    rgb(:,:,3) = 0;
    figure(fig),imshow(rgb,[]);
end