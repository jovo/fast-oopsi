function result = CTassignLinkColor(clink, ccolor)
[ti tj] = find(clink ~= 0);
uti = unique(ti);
count = 0;
result = [];
for i = 1 : length(uti)
    if isempty(ccolor)
        c = uint8(rand(1,3)*255);
    else
        iindex = find(ccolor{1}(:,1) == uti(i));
        if isempty(iindex)
            c = uint8(rand(1,3)*255);
        else
            c = ccolor{1}(iindex(1),2:4);
        end
    end
    result{1}(i,:) = [uti(i) c];
    jindex = find(ti == uti(i));
    for j = 1: length(jindex)
        count = count + 1;
        result{2}(count,:) = [tj(jindex(j)) c];
    end
end