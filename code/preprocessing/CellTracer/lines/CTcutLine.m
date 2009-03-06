function result = CTcutLine(bw,MergeThreshold, SplitThreshold, MergeLineSize, SplitLineSize, FilterLineSize, maxcellwidth)
ballsize = 1;
notdone = 1;
while notdone>0
    temp = imopen(bw,strel('disk',ballsize));
    %figure(101),imshow(temp,[]);pause
    [L num] = bwlabel(temp,4);
    if num ==1
        ballsize = ballsize + 1;
    else
        notdone = 0;
        ballsize = ballsize - 1;
    end
end
if ballsize > 0
    im = imopen(bw,strel('disk',ballsize));
else
    im = bw;
end

center = bwmorph(im,'thin',Inf); 
figure(101),imshow(im);
figure(105),imshow(center);pause
[lines splitmask] = CTgetLineSegments(center);

MergeSplit = 1; maxiter= 100;
while MergeSplit && maxiter > 0;
    %maxiter
    [lines scores split] = CTsplitLines(lines,SplitLineSize, SplitThreshold);
    [lines merge] = CTmergeLines(lines, splitmask,MergeLineSize, MergeThreshold);
    MergeSplit = split >0 | merge > 0;
    maxiter = maxiter - 1;
end

[m n] = size(bw);
newlines = CTfilterLines(m,n,lines,FilterLineSize);
lines = newlines;
center = CTreconstructCenterFromLines(m,n,lines); 
area = bwmorph(center,'thicken',maxcellwidth);
bw(area ==0) = 0;
result = bw;
figure(102),imshow(result,[]);