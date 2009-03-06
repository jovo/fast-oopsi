function result = CTExpand(im,mask,conn,refdisk,greedy)
%h1 = figure;h2 = figure;
if isempty(mask)
    result  = im ; 
else
    result = bwlabel(mask);
end
if ~isempty(refdisk) && refdisk > 0
    temp = 1-mask; temp = imdilate(temp,strel('disk',refdisk));
    im(temp==0) = 0;
    temp = im(mask>0 & temp >0); 
else
    temp = im(mask>0);
end
vmin = min(temp(:)); vmax = max(temp(:));

if ~isempty(greedy) && greedy > 0
    temp = imdilate(mask,strel('disk',greedy));
    im(mask == 0 & temp ==0) = 255;
end
for i=vmin:vmax 
    b = im<=i;    
    [L num] = bwlabel(b,conn);
    for k = 1:num
        currentobject = result; currentobject(L ~= k) = 0;
        t1  = unique(currentobject(:)); tmin = min(t1); tmax = max(t1);
        if length(t1) == 2 && tmin==0  %expand the object
            result(L==k) = tmax;
        elseif length(t1) == 2 && tmin > 0 %merge without add new points
            %result (L==k) = tmin;
        else %objects start to merge together
            notdone = length(currentobject(currentobject>0));
            while notdone
                for m=1:length(t1)
                    if t1(m)>0
                        marker = result < 0;
                        marker(result == t1(m)) = 1;
                        marker = imdilate(marker,strel('disk',1));
                        result(marker > 0 & L ==k ) = t1(m);
                        currentobject(marker > 0 & L ==k) = 0;
                    end
                end
                notdone = length(currentobject(currentobject>0));
            end
        end
    end
    rgb(:,:,1) = im;rgb(:,:,2) = uint8(mask)*255;rgb(:,:,3) = uint8(result) * 255;
    figure(100),imshow(rgb)
end