function conn = CTgetConnections(im)
[m n] = size(im);
conn = zeros(m,n,'uint8');
for i = 2: m-1
    for j= 2:n-1
        if im(i,j) > 0
            temp = im(i-1:i+1,j-1:j+1);
            conn(i,j) = sum(temp(:)) - 1;
        end
    end
end