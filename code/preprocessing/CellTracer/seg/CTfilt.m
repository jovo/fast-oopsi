function [Mean,newIm] = CTfilt(im, ballsize,method,standardize,percentile, hint, hint2)
[m n] = size(im);
Mean = zeros(m,n);

se = strel('disk',ballsize+1);
ball= getnhood(se);

ballarea = double(sum(ball(:)));
if standardize == 0
    Mean = double(im);
elseif standardize == 1
     %estat = erangefilt(im,ballsize);
     estat = CTblock(im,[],[],ballsize, {@mean},[]);
     Mean = double(im) - estat(:,:,1);
     Mean = Mean - min(Mean(:));
     Mean(im == 0) = 0;
elseif standardize == 2
    estat = CTblock(im,[],[],ballsize, {@min},[]);
    Mean = double(im) - estat(:,:,1);
    Mean = Mean - min(Mean(:));
    Mean(im == 0) = 0;
end
if method == -1
    newIm = Mean;
    CutThreshold = percentile;
    newIm(newIm >= CutThreshold) = 255;
    newIm(newIm<CutThreshold) = 0;
    newIm = uint8(newIm);
elseif method == 0
    newIm = Mean;
    temp1 = newIm(im>0);
    CutThreshold = ceil(prctile(temp1,percentile));
    newIm(newIm >= CutThreshold) = 255;
    newIm(newIm<CutThreshold) = 0;
    newIm = uint8(newIm);
elseif method == 1
    newIm = CTblock(Mean,[],[],ballsize,{@method1},{percentile,ballarea,hint});
    if hint2 <= 0
        newIm(newIm>0) = 255;
    else
        temp = double(newIm(newIm>0)); pert = ceil(prctile(temp,hint2));
        newIm(newIm>= pert) = 255;
        newIm(newIm < pert) = 0;
    end
    newIm = uint8(newIm);
elseif method == 2
    newIm = CTblock(Mean,[],[],ballsize,{@method2},{percentile});
    if hint2 <= 0
        newIm(newIm>0) = 255;
    else
        temp = double(newIm(newIm>0)); pert = ceil(prctile(temp,hint2));
        newIm(newIm>= pert) = 255;
        newIm(newIm < pert) = 0;
    end
    newIm = uint8(newIm);
elseif method == 3
    %get a bw image first
    border =im;border(border>0) = 255; border = bwperim(border);
    %figure(100),imshow(Mean,[]);
    newIm = CTblock(Mean,[],border,ballsize,{@method3},{percentile,ballarea,hint});
    if hint2 <= 0
        newIm(newIm>0) = 255;
    else
        temp = newIm; 
        temp = temp(temp>0); 
        if ~isempty(temp)
            pert = ceil(prctile(temp,hint2));
            newIm(newIm>= pert) = 255;
            newIm(newIm < pert) = 0;
        end
    end
    newIm = uint8(newIm);
elseif method == 4
    %get a bw image first
    border =im;border(border>0) = 255;
    border = bwperim(border,8);
    newIm = border;
    newIm = uint8(newIm);
else
    %not defined yet
end
