function newim = CTreconstructCenterFromLines(m,n, lines)
newim = zeros(m,n,'uint8');
for j = 1:length(lines)
    currentline = lines{j}.points;
    [pm pn] = size(currentline);
    for k=1:pm
        newim(currentline(k,1), currentline(k,2)) = 1;
    end
end