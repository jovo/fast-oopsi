function result = CTblock(im,inputmask,outputmask,ballsize,fun,para)
[m n] = size(im);
if isempty(inputmask) 
    inputmask = im; inputmask(inputmask > 0) = 255;
end
if isempty(outputmask) 
    outputmask = im; outputmask(outputmask>0) = 255;
end
if isempty(fun)
    fun{1} = @max;
    fun{2} = @min;
    fun{3} = @mean;
end
im(inputmask == 0) = 0;

impadding = zeros(m + 2 * ballsize, n + 2 * ballsize,'uint8');
impadding(ballsize+1:ballsize+m,ballsize+1:ballsize+n) = im;

if ballsize<3
    se = strel('disk',ballsize);
else
    se = strel('disk',ballsize+1); 
end
ball= getnhood(se);
ball = ball > 0;

result = zeros(m,n,length(fun));
for i=ballsize+1:ballsize+m
    indexi = i-ballsize;
    for j=ballsize+1:ballsize+n
        indexj = j-ballsize;
        if outputmask(indexi,indexj) > 0
            temp1 = impadding(i-ballsize:i+ballsize, j-ballsize:j+ballsize);
            temp1 = temp1(ball & temp1 >0); 
            if ~isempty(temp1)
                for k = 1:length(fun)
                    if isempty(para)
                        result(indexi,indexj,k) = fun{k}(double(temp1));
                    else
                        result(indexi,indexj,k) = fun{k}(temp1,impadding(i, j),para);
                    end
                end
            end
        end
    end
end
