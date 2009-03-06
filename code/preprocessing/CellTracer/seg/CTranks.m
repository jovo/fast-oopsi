function result = CTranks(x,v,para)
if para{1} == 1
    result = length(find(x<v)) / length(x);
elseif para{1} == 2
    result = length(find(x>v)) / length(x);
elseif para{1} == 3
    y= unique(x);
    result = length(find(y<v)) / length(y);
elseif para{1} == 4
    y= unique(x);
    result = length(find(y>v)) / length(y);
elseif para{1} == 5
    xmin = min(x); xmax = max(x);
    if xmin == xmax
        result = 0;
    else
        result = double(v-xmin)/double(xmax - xmin);
    end
elseif para{1} ==6
    xmin = min(x); xmax = max(x);
    if xmin == xmax
        result = 1;
    else
        result = double(xmax-v)/double(xmax - xmin);
    end
elseif para{1} ==7 
    result = length(find(x==v)) / length(x);
elseif para{1} == 8
    y= unique(x);
    result = length(find(y==v)) / length(y);
end
