function [newlines,scores split] = CTsplitLines(lines, linesize, MergeSplitThreshold)
newlines = {}; nCount = 1; scores = [];
split = 0;
halfsize = floor(linesize/2);
if halfsize * 2 + 1 > linesize 
    halfsize = halfsize - 1;
end

for k=1:length(lines)
    if length(lines{k}.points) == 0
        continue;
    end
    
    [m n] = size(lines{k}.points);
    if m < linesize
        newlines{nCount}.points = lines{k}.points;
        nCount = nCount + 1;
        continue;
    end
    scores = [];
    for i = halfsize+1:m-halfsize-1
        seg.points = lines{k}.points(i-halfsize:i+halfsize,:);
        score = CTgetLineFitScore(seg.points);
        scores = [scores; [i score]];
    end
    
    maxscore = max(scores(:,2));
    if maxscore > MergeSplitThreshold
        SplitIndex = find(scores(:,2) == maxscore);
        SplitIndex = scores(SplitIndex,1);
        newlines{nCount}.points = lines{k}.points(1:SplitIndex,:);
        nCount = nCount + 1;
        newlines{nCount}.points = lines{k}.points(SplitIndex:end,:);
        nCount = nCount + 1;
        split = 1;
    else 
        newlines{nCount}.points = lines{k}.points;
        nCount = nCount + 1;
    end
end

    
