function result = CTshiftImage(hObject,mask,frame)
%get the alignment
[clink ccolor alignment] = getCurrentLink(hObject,frame);
if isempty(alignment)
    result = mask;
else
    [m n] = size(alignment.score); shift = (m+1)/2;
    [m n] = size(mask);
    bestlocation = alignment.bestlocation;
    container = zeros(m+2* shift, n+ 2 * shift,'uint8');
    t2 = container; t2(bestlocation(1):bestlocation(1) + m-1, bestlocation(2):bestlocation(2)+n-1) = mask;
    result = t2(shift+1:shift+m, shift+1:shift+n);
    result = result > 0;
    %figure(100),imshow(result);
end