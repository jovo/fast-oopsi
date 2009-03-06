function [newline,score] = CTgetLineMergeScore(line1, line2, focuspoints,linesize)
[m n] = size(focuspoints);
halfsize = floor(linesize/2);
if n == 2
    if halfsize * 2 + 1 > linesize %just make sure we have enough points
        halfsize = halfsize - 1;
    end
    point1 = focuspoints;
    point2 = point1;
    newline= point1;
else
    if halfsize * 2 + 2 > linesize %just make sure we have enough points
        halfsize = halfsize - 1;
    end
    point1 = focuspoints(1,1:2);
    if n==4
        point2 = focuspoints(1,3:4);
        newline= [point1; point2];
    else
        point2 = focuspoints(1,5:6);
        newline= [point1; focuspoints(1,3:4); point2];
    end
end

%figure out the points used from each line.
size1 = length(line1(:,1)) -1; size2 = length(line2(:,1))-1;
if size1 < halfsize
    size2 = halfsize * 2 - size1;
elseif size2 < halfsize
    size1 = halfsize * 2 - size2;
else
    size1 = halfsize;
    size2 = halfsize;
end

fullline = newline;
if point1(1,1) == line1(1,1) && point1(1,2) == line1(1,2)
    newline = [newline;line1(2:size1+1,:)];
    fullline = [fullline; line1(2:end,:)];
else
    [m n] = size(line1);   %something wrong here, so a lazy fix for now
    for i = 1:min(size1, m-1);
        newline = [newline;line1(m-i,:)];
    end
    for i = 1:length(line1(:,1)) - 1;
        fullline = [fullline;line1(m-i,:)];
    end
end
if point2(1,1) == line2(1,1) && point2(1,2) == line2(1,2)
    [m n] = size(line1); %something wrong here, so a lazy fix for now
    for i=1:min(size2,m-1);
        newline = [line2(1+i,:);newline];
    end
    for i=1:1:length(line2(:,1)) - 1;
        fullline = [line2(1+i,:);fullline];
    end
else
    [m n] = size(line2);
    if m>size2 && m > 1 && size2 >= 1 
        newline = [line2(m-size2:end-1,:);newline];
        fullline = [line2(1:m-1,:);fullline];
    end
    
end

score = CTgetLineFitScore(newline);
newline = fullline;
