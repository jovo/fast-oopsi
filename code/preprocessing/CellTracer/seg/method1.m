function result = method1(x,v,para)
result = 0;
localratio = double(para{2} - length(x)) / para{2};
if localratio < para{3}
    %we modify the percentile here a bit in this case
    %the formula is kind of for no reason here
    localpercentile = max(para{1} - (para{3} - localratio) / 2, para{1} * 0.75);
    CutThreshold = ceil(prctile(x,localpercentile));
    if v >= CutThreshold 
         result = v;
    end
end
