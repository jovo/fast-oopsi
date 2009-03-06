function result = method3(x,v,para)
result = 0;
localratio = double(para{2} - length(x)) / para{2};
if localratio < para{3}
   CutThreshold = ceil(prctile(x,para{1}));
    if v > CutThreshold
        result = v;
    end
end
                