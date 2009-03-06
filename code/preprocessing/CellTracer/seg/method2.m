function result = method2(x,v,para)
result = 0;
if length(x) >= 5 
    CutThreshold = ceil(prctile(x,para{1}));
    if v >= CutThreshold 
         result = v;
    end
end