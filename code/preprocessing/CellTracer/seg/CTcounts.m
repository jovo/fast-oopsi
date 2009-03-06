function result = CTcounts(x,v,para)
t1 = length(find(x<para{1}));
t2 = length(find(x>para{2}));
if t2 ==0 
    result = -1;
else
    result = t1 / t2;
end
result(result == -1) = max(result(:));



