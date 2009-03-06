function alignment = CTglobalAlignment(im1,im2,para,debug)
shift = para{1};
%evaluate all possible shiftings
[m n] = size(im1);
container = zeros(m+2* shift, n+ 2 * shift);
resultsize = 2 * shift+1;
result = zeros(resultsize, resultsize);
locations = zeros(resultsize, resultsize);
t1 = container; t1(shift+1:shift+m, shift+1:shift+n) = im1;
for i = 1:resultsize
    for j = 1:resultsize
        t2 = container; t2(i:i + m-1, j:j+n-1) = im2;
        bw = double(t1>0 & t2>0); 
        result(i,j) = sum(double(bw(:)));      
        locations(i,j) = double((i-shift) * (i-shift) + (j-shift) * (j-shift));
    end
end
locations = locations+0.1; % all locations should be non-zero
%find the best location
bestshift = result; bestshift(bestshift < max(bestshift(:))) = 0;
locations(bestshift==0) = 0;
if length(find(locations>0))>1 %if more than one best locations found 
    temp = locations(locations>0);
    if length(find(locations>min(temp(:)))) < length(find(locations>0))
        locations(locations>min(temp(:))) = 0;
    end
end
alignment.score = result;
try 
    [r c] = find(locations>0); bestlocation = [r(1) c(1)];
catch
    bestlocation = [shift+1 shift+1];
end
alignment.bestlocation = bestlocation;
alignment.relativeshift = bestlocation-[shift+1 shift+1];
if debug > 0
    CTshowGlobalAlignment(im1,im2,alignment,100);
end