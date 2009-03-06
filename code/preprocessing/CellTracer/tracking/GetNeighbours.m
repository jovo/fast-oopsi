%this fuction is very much like BuildLinks,but used to find the neighbours
%for each cell
function result = GetNeighbours(label,d1,d2,method)
[m n] = size(label);
n1 = max(label(:));
result = zeros(n1,n1);
d = d1 * d2;
if d>0
    if method == 1 %to be fixed later to allow deleted cells!!!!!!
        container = zeros(m+2* d, n+ 2 * d,'uint8');
        t1 = container; t1(d+1:d+m, d+1:d+n) = label; 
        for i = 1:2 * d1+1
            row = (i-1) * d2 + 1;
            for j = 1:2* d1+1
                col = (j - 1 ) * d2 + 1;
                t2 = container; t2(row:row + m-1, col:col+n-1) = label;
                bw = t1 & t2;
                [L num] = bwlabel(bw,4);

                num1 = zeros(n1,n1);

                for j = 1:num
                    [r c] =find(L == j);
                    commpix = length(r);
                    index1= t1(r(1),c(1));
                    index2 =t2(r(1),c(1)); 
                    num1(index1, index2) = num1(index1, index2) + commpix;
                end
                result = result + num1; 
            end
        end
        result(result>0) = 1;
    elseif method == 2 %to be fixed later to allow deleted cells!!!!!!
        boundingbox = zeros(n1,4);
        for i = 1:n1
            [r c] = find(label == i);
            boundingbox(i,1) = max(1,min(r)-d);
            boundingbox(i,2) = max(1,min(c)-d);
            boundingbox(i,3) = min(m,max(r)+d);
            boundingbox(i,4) = min(n,max(c)+d);
        end
        for i = 1:n1
            for j=1:i
                if boundingbox(j,1) >= boundingbox(i,1) && boundingbox(j,1) <= boundingbox(i,3) && boundingbox(j,2) >= boundingbox(i,2) && boundingbox(j,2) <= boundingbox(i,4)
                    result(i,j) = 1; result(j,i) = 1;
                elseif boundingbox(j,1) >= boundingbox(i,1) && boundingbox(j,1) <= boundingbox(i,3) && boundingbox(j,4) >= boundingbox(i,2) && boundingbox(j,4) <= boundingbox(i,4)
                    result(i,j) = 1; result(j,i) = 1;
                elseif boundingbox(j,3) >= boundingbox(i,1) && boundingbox(j,3) <= boundingbox(i,3) && boundingbox(j,2) >= boundingbox(i,2) && boundingbox(j,2) <= boundingbox(i,4)
                    result(i,j) = 1; result(j,i) = 1;
                elseif boundingbox(j,3) >= boundingbox(i,1) && boundingbox(j,3) <= boundingbox(i,3) && boundingbox(j,4) >= boundingbox(i,2) && boundingbox(j,4) <= boundingbox(i,4)
                    result(i,j) = 1; result(j,i) = 1;
                end 
            end
        end
    else   %method == 3
        boundingbox = zeros(n1,4);
        for i = 1:n1
            [r c] = find(label == i);
            if ~isempty(r)
                boundingbox(i,1) = max(1,min(r)-d);
                boundingbox(i,2) = max(1,min(c)-d);
                boundingbox(i,3) = min(m,max(r)+d);
                boundingbox(i,4) = min(n,max(c)+d);
            end
        end
        for i = 1:n1
            if any(boundingbox(i,:))
                t2 = label(boundingbox(i,1):boundingbox(i,3),boundingbox(i,2):boundingbox(i,4));
                t1 = t2 * 0;
                t1(t2 == i) = 255;
                t1 = uint8(bwmorph(t1,'thicken',d));
                bw = t1 & t2;
                [L num] = bwlabel(bw,4);
                for j = 1:num
                    [r c] =find(L == j);
                    index =t2(r(1),c(1)); 
                    result(i,index) = 1;
                    result(index,i) = 1;
                end
            end
        end
    end
else
    for i = 1:n1
        result(i,i) = 1;
    end
end
   
