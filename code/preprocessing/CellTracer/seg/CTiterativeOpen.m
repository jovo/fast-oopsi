function result = CTiterativeOpen(im,disk,recursive)
im = im>0;
ime = imopen(im,strel('disk',disk));
result = im; result(ime==0) = 0;
if recursive>0 
    Diff= im - result;
    if any(Diff>0)
        result = CTiterativeOpen(result,disk,recursive);
    end
end    

