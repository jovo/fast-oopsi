function [newlines merge]= CTmergeLines(lines, SplitMask, linesize, MergeThreshold)
%merge the lines
numberoflines = length(lines);
splitpoints = find(SplitMask>0);
if length(lines) > 0
    canmerge = 1;
    while canmerge > 0
        canmerge = 0;        
        %first merge lines with common ends
        %line a, b to be merged only if they have common ends and have a
        %combined length of at least linesize long
        
        %for each line, find lines that have common ends and combined size
        %of at least linesize
        for i=1:length(lines)
            if length(lines{i}.points) > 0
                mergecandidate = [];
                commonpoints = [];
                currentlinelength = length(lines{i}.points);
                pt1 = lines{i}.points(1,:);
                pt2 = lines{i}.points(end,:);
                for j = 1:length(lines)
                    if j ~= i && length(lines{j}.points) > 0
                        if length(lines{j}.points) + currentlinelength > linesize
                            pt3 = lines{j}.points(1,:);
                            pt4 = lines{j}.points(end,:);
                            if pt1(1) == pt3(1) && pt1(2) == pt3(2)
                                mergecandidate = [mergecandidate j];
                                commonpoints = [commonpoints; pt1];
                            elseif pt1(1) == pt4(1) && pt1(2) == pt4(2) 
                                mergecandidate = [mergecandidate j];
                                commonpoints = [commonpoints; pt1];
                            elseif pt2(1) == pt3(1) && pt2(2) == pt3(2)
                                mergecandidate = [mergecandidate j];
                                commonpoints = [commonpoints; pt2];
                            elseif pt2(1) == pt4(1) && pt2(2) == pt4(2)
                                mergecandidate = [mergecandidate j];
                                commonpoints = [commonpoints; pt2];
                            end
                        end
                    end
                end
                if length(mergecandidate) > 0 %we have some possible matchs here
                    %find the best candidate
                    BestScore = 100; BestLine = 0; BestCandidate = 0; 
                    for j = 1:length(mergecandidate)
                        [newline score] =  CTgetLineMergeScore(lines{i}.points, lines{mergecandidate(j)}.points, commonpoints(j,:), linesize);
                        if score < BestScore
                            BestScore = score; BestLine = newline; BestCandidate = j; 
                        end
                    end
                    if  BestScore < MergeThreshold && BestCandidate > 0
                        lines{i}.points = BestLine;
                        lines{mergecandidate(BestCandidate)}.points = [];
                        canmerge = 1;  %so we made some merges here and need to rerun until no more merge can be made
                    end
                end
            end
        end
        
        %now deal with merger of lines close touch but have no common ends
        for i=1:length(lines)
            if length(lines{i}.points) > 0
                mergecandidate = [];
                touchpoints = [];
                currentlinelength = length(lines{i}.points);
                pt1 = lines{i}.points(1,:);
                pt2 = lines{i}.points(end,:);
                for j = 1:length(lines)
                    if j ~= i && length(lines{j}.points) > 0
                        if length(lines{j}.points) + currentlinelength >= linesize
                            pt3 = lines{j}.points(1,:);
                            pt4 = lines{j}.points(end,:);
                            if all(abs(pt1 - pt3) <=1) && any(abs(pt1 - pt3)> 0) 
                                mergecandidate = [mergecandidate j];
                                touchpoints = [touchpoints; [pt1 pt3]];
                            elseif all(abs(pt1 - pt4) <=1) && any(abs(pt1 - pt4)> 0)  
                                mergecandidate = [mergecandidate j];
                                touchpoints = [touchpoints; [pt1 pt4]];
                            elseif all(abs(pt2 - pt3) <=1) && any(abs(pt2 - pt3)> 0)  
                                mergecandidate = [mergecandidate j];
                                touchpoints = [touchpoints; [pt2 pt3]];
                            elseif all(abs(pt2 - pt4) <=1) && any(abs(pt2 - pt4)> 0)  
                                mergecandidate = [mergecandidate j];
                                touchpoints = [touchpoints; [pt2 pt4]];
                            end
                        end
                    end
                end
                if length(mergecandidate) > 0 %we have some possible matchs here
                    %find the best candidate
                    BestScore = 100; BestLine = 0; BestCandidate = 0; 
                    for j = 1:length(mergecandidate)
                        [newline score] =  CTgetLineMergeScore(lines{i}.points, lines{mergecandidate(j)}.points, touchpoints(j,:), linesize);
                        if score < BestScore
                            BestScore = score; BestLine = newline; BestCandidate = j; 
                        end
                    end
                    if BestScore < MergeThreshold
                        lines{i}.points = BestLine;
                        lines{mergecandidate(BestCandidate)}.points = [];
                        canmerge = 1;  %so we made some merges here and need to rerun until no more merge can be made
                    end
                end
            end
        end
        
        %now deal with merger of lines no touching but can be connected by
        %a split point
        for i=1:length(lines)
            if length(lines{i}.points) > 0
                mergecandidate = [];
                tripoints = {}; nCount = 1;
                currentlinelength = length(lines{i}.points);
                pt1 = lines{i}.points(1,:);
                pt2 = lines{i}.points(end,:);
                for j = 1:length(lines)
                    if j ~= i && length(lines{j}.points) > 0
                        if length(lines{j}.points) + currentlinelength >= linesize
                            pt3 = lines{j}.points(1,:);
                            pt4 = lines{j}.points(end,:);
                            pts = CTgetBridgePoint(pt1,pt3,splitpoints);
                            if length(pts) > 0
                                mergecandidate = [mergecandidate j];
                                tripoints{nCount} = pts;
                                nCount = nCount+1;
                            else
                                pts = CTgetBridgePoint(pt1,pt4,splitpoints);
                                if length(pts) > 0
                                    mergecandidate = [mergecandidate j];
                                    tripoints{nCount} = pts;
                                    nCount = nCount+1;
                                else
                                    pts = CTgetBridgePoint(pt2,pt3,splitpoints);
                                    if length(pts) > 0 
                                        mergecandidate = [mergecandidate j];
                                        tripoints{nCount} = pts;
                                        nCount = nCount+1;
                                    else
                                        pts = CTgetBridgePoint(pt2,pt4,splitpoints);
                                        if length(pts) > 0
                                            mergecandidate = [mergecandidate j];
                                            tripoints{nCount} = pts;
                                            nCount = nCount+1;
                                        end
                                    end
                                end
                            end    
                        end
                    end
                end
                if length(mergecandidate) > 0 %we have some possible matchs here
                    %find the best candidate
                    BestScore = 10; BestLine = 0; BestCandidate = [0 0]; 
                    for j = 1:length(mergecandidate)
                        pts = tripoints{j};
                        for k = 1:length(tripoints)
                            [newline score] =  CTgetLineMergeScore(lines{i}.points, lines{mergecandidate(j)}.points, pts{k}, linesize);
                            if score < BestScore
                                BestScore = score; BestLine = newline; BestCandidate = [j k]; 
                            end    
                        end
                        
                    end
                    if BestScore < MergeThreshold
                        lines{i}.points = BestLine;
                        lines{mergecandidate(BestCandidate(1))}.points = [];
                        canmerge = 1;  %so we made some merges here and need to rerun until no more merge can be made
                    end
                end
            end
        end
    end
end
newlines = {}; count = 0;
for i = 1:length(lines)
    if length(lines{i}.points) > 0;
        count = count + 1;
        newlines{count} = lines{i};
    end
end
if numberoflines ~= length(newlines)
    merge = 1;
else
    merge = 0;
end