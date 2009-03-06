function points = CTgetBridgePoint(pt1, pt2, splitlinks)
[m n] = size(splitlinks);
points = {}; nCount = 1;
if all(abs(pt1 - pt2) <=2) && any(abs(pt1 - pt2)> 1)
    for i=1:m
        pt = splitlinks(i,:);
        if all(abs(pt1 - pt) <=1) && any(abs(pt1 - pt)> 0)
            if all(abs(pt2 - pt) <=1) && any(abs(pt2 - pt)> 0)
                points{nCount} = [pt1 pt pt2];
                nCount = nCount + 1;
            end
        end
    end
end
