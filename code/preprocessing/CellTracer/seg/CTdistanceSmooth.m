function result = CTdistanceSmooth(im,dist,tol,debug)
bw = im==0;

%first do a distance transform
d1 = bwdist(bw);
if dist ==0
    d = d1;
else
    bw = d1<=dist;
    d = bwdist(bw);
end
% if debug > 0
%     figure(100),imshow(d1,[]);
%     figure(101),imshow(d,[]);
% end

a = d>0;
if tol > 0
    if tol < 1
        tol = tol * 100;
    end
    [L num] = bwlabel(a);
    for i=1:num
        temp = d(L==i);
        if length(temp) > 5
            CutThreshold = prctile(temp,tol);
            if length(find(temp<=CutThreshold)) / length(temp) * 100 <= tol + 5
                d(L==i & d<=max(CutThreshold,1)) = 0;
            end
        else
            d(L == i) = 0;
        end
    end
else
    d = a;
end
im(d==0) = 0;
result = im;