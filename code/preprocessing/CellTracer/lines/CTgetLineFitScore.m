function score = CTgetLineFitScore(line)
x = line(:,1); y = line(:,2);
n = length(x);
x1n = x(n)-x(1); y1n = y(n)-y(1);
dd = sqrt(x1n*x1n + y1n*y1n);

%sum = 0;
%for i=2:n-1
%    sum = sum + abs(x1n * (y(1)-y(i)) - (x(1)-x(i)) * y1n);
%end
%score = sum / (dd * n);
%Another way
max = 0;
for i=2:n-1
    temp = abs(x1n * (y(1)-y(i)) - (x(1)-x(i)) * y1n);
    if temp > max
        max = temp;
    end
end
score = max / dd ;
