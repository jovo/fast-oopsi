function result = ExpandFrame(a)
a = double(a);
[m n] = size(a);
b = zeros(2 * m -1, 2 * n -1);
for i=1:m
    for j = 1:n
        b(2 * i -1, 2 * j -1) = a(i,j);
    end
end
for i=1:m-1
    for j = 1:n-1
        b(2 * i -1, 2 * j) = (a(i,j) + a(i,j+1))/2;
        b(2 * i, 2* j -1) = (a(i,j) + a(i+1,j))/2;
        b(2*i, 2*j) = (a(i,j)+a(i,j+1)+ a(i+1,j) + a(i+1,j+1))/4;
    end
end

for i=1:n-1
    b(2*m-1,2*i) = (a(m,i) + a(m,i+1))/ 2;
end

for i=1:m-1
    b(2*i,2*n-1) = (a(i,n) + a(i+1,n))/ 2;
end
result = uint8(b);