function newlines = CTfilterLines(m,n,lines, linesize) 
%First get ride of short lines
HalfLineSize = floor(linesize / 2);
HalfHalfLineSize = floor(HalfLineSize / 2);
MaxLineWidth = HalfLineSize;

%remove short lines
newlines = {}; count = 0; 
for i = 1:length(lines)
    if length(lines{i}.points) >= HalfHalfLineSize;
        count = count + 1;
        newlines{count} = lines{i};
    end
end
lines = newlines;

%Trim lines that are connected to other lines
CenterIm = CTreconstructCenterFromLines(m,n, lines);
SplitMask = CTgetConnections(CenterIm);
for i = 1:length(lines)
    [m n] = size(lines{i}.points);
    x = lines{i}.points(1,1);
    y = lines{i}.points(1,2);
    if SplitMask(x,y) == 2
        lines{i}.points = lines{i}.points(min(HalfHalfLineSize + 1,m):m,:);
    elseif SplitMask(x,y) > 2
        lines{i}.points = lines{i}.points(min(HalfLineSize + 1,m):end,:);
    end
    
    [m n] = size(lines{i}.points);
    x = lines{i}.points(m,1);
    y = lines{i}.points(m,2);
    if SplitMask(x,y) == 2
        lines{i}.points = lines{i}.points(1:max(1,m-HalfHalfLineSize),:);
    elseif SplitMask(x,y) > 2
        lines{i}.points = lines{i}.points(1:max(1,m-HalfLineSize),:);
    end
end

%remove short lines again
newlines = {}; count = 0;
for i = 1:length(lines)
    if length(lines{i}.points) >= HalfLineSize;
        count = count + 1;
        newlines{count} = lines{i};
    end
end
