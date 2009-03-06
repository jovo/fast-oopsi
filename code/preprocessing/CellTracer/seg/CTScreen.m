function [result result2] = CTScreen(a,mask,maxsize,minsize,disk,conn,limit,maxspan,hint,para)
%h1 = figure;h2 = figure;
abackup = a;
if isempty(mask)
    result  = zeros(size(a)); 
else
    [result count] = bwlabel(mask);
end
result2 = result;
if isempty(limit)
    vmin = min(a(:)); vmax = max(a(:));
    venter = vmin; vexit = vmax;
else
    vmin = limit(1); venter = limit(2); vexit = limit(3); vmax = limit(4);
    if vmin < 0 
        vmin = min(a(:));
    end
    if venter < 0 
        venter = vmin;
    end
    if vmax < 0
        vmax = max(a(:)); 
    end
    if vexit < 0
        vexit = vmax;
    end
    a(a<vmin) = vmin; a(a>vmax) = vmax;
end 
indexes = sort(unique(a(:)));
dcount = 0;
for i = 1:length(indexes) 
    b = a<indexes(i);
    %figure(22),imshow(b);
    
    %size and intensity value constraints
    c = CTiterativeApply(double(b),0, conn,{@(x) double(x>0).* double(length(x(x>0)))}, []);
     c(c<minsize) = 0;
    if isempty(maxspan)
        c(c>maxsize) = 0;
    else
        c(c>maxsize*2) = 0;
    end
    
    %refine result2
    temp = c>0 & result2 > 0; temp = imreconstruct(temp,c>0);
    result2(temp>0) = 255;
    
    %smoothing if necessary
    if disk >0
        d  = c>0;d1 = imopen(d,strel('disk',disk));
        d1 = CTiterativeApply(double(d1),0, conn,{@(x) double(x>0).* double(length(x(x>0)))}, []);
        d1(d1<minsize) = 0; d1 = d1>0;
        d = imreconstruct(d1,d,conn);
        if disk-1>0
            d = imopen(d,strel('disk',disk-1));
            d = CTiterativeApply(double(d),0, conn,{@(x) double(x>0).* double(length(x(x>0)))}, []);
            d(d<minsize) = 0; d = d>0;
        end
    else
        d = c>0;
    end
    
    [L num] = bwlabel(d,conn);
    if num > 0
        for k = 1:num
            currentobject = result; currentobject(L ~= k) = 0;
            t1  = unique(currentobject(:)); tmin = min(t1);tmax = max(t1);
            if length(t1) == 1 %first time object
                if isempty(maxspan) || length(find(L==k)) < maxsize
                    t2 = a(L==k); t2min = min(t2(:)); t2max = max(t2(:));
                    if t2min >= venter && t2max <= vexit
                        count = count + 1;
                        result(L==k) = count;
                        objspan(count) = i;
                        result2(L==k) = 255;
                    end
                end
            elseif length(t1) == 2
                if tmin == 0   %expand the object
                    if ~isempty(maxspan) % if there is maxspan contraint, expand the way defined below
                        if i >objspan(max(t1)) + maxspan && disk > 0
                            temp = L==k; 
                            temp(result == tmax) = 0;
                            d1 = imopen(temp,strel('disk',disk));
                            d1 = CTiterativeApply(double(d1),0, conn,{@(x) double(x>0).* double(length(x(x>0)))}, []);
                            d1(d1<minsize) = 0; d1 = d1>0;
                            [L1 num1] = bwlabel(d1,conn);
                            if num1 > 0
                                for m =1:num1
                                    if length(find(L1==m)) < maxsize
                                        count = count+1;
                                        result(L1==m) = count;
                                        objspan(count) = i;
                                    end
                                end
                            end
                        end
                    else % just expand the object
                        result(L==k) = tmax;
                    end
                else %can this happen?
                    msgbox('strange case, please report to the author');
                end
            else  %objects start to merge together
                currentobject = result; currentobject(L ~= k) = 0;
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
        if (dcount==10)
            rgb = label2rgb(result,'spring','c','shuffle');
            figure(100),imshow(rgb);
            dcount = 0;
        else
            dcount = dcount + 1;
        end
    end
end
if ~isempty(hint) 
    if hint == 1 
        packratio = para;
        for i = 1:count
            temp = result == i;
            a = double(sum(temp(:)));
            b = bwperim(temp,conn); b = double(sum(b(:)));
            if b / a > packratio
                result(result == i) = 0;
            end
        end
    elseif hint == 2 %mean background
        scale = para;
        for i = 1:count
            temp = result == i;
            [x y] = find(temp>0);
            bb = abackup(min(x):max(x),min(y):max(y));
            btemp = temp(min(x):max(x),min(y):max(y));
            a = bb(btemp==0); amean = mean(double(a)); astd = sqrt(var(double(a)));
            b = bb(btemp>0);  bmean = mean(double(b));
            if bmean>=amean-scale * astd && bmean<= amean+scale * astd
                result(result == i) = 0;
            end
        end
    elseif hint ==3 %min value
        minvalue = para;
        for i = 1:count
            temp = result == i;
            temp = abackup(temp>0); 
            if min(temp(:)) > minvalue
                result(result == i) = 0;
            end
        end
    elseif hint == 4 %min value,percentile
        minvalue = vexit;
        ratio = para;
        for i = 1:count
            temp = abackup(result == i); 
            cratio = double(length(temp(temp<=minvalue))) / double(length(temp(temp>0)));
            if cratio < ratio
                result(result == i) = 0;
            end
        end
    end
    %rgb = label2rgb(result,'spring','c','shuffle');
    %figure(100),imshow(rgb);
    
end