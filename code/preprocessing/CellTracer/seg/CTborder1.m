function bordermask = CTborder1(im,mask,para,debug)
cminmax = ones(size(im));
cmaxmin = zeros(size(im));
for ballsize = 1:para{4}
    se = strel('disk',ballsize);
    erange(:,:,1) = imdilate(im,se);
    erange(:,:,2) = imerode(im,se);
    temp = double(erange(:,:,1) - erange(:,:,2)); 
    
    temp2 = double(im - erange(:,:,2));
    temp2(temp>0) = temp2(temp>0) ./ temp(temp>0);
    cminmax = min(cminmax,temp2);
    
     temp2 = double(erange(:,:,1)-im);
     temp2(temp>0) = temp2(temp>0) ./ temp(temp>0);
     cmaxmin = max(temp2,cmaxmin);
end

%deal with the flat regions
temp = cmaxmin + cminmax;
temp1 = cminmax == 1 | temp == 0;
temp2 = temp ==0;
temp = imreconstruct(temp2,temp1,8);
cminmax(temp>0) = 1.0;

if debug > 0
    figure(100),imshow(cminmax,[]);
    %figure(101),imshow(cmaxmin,[]);
end

bordermask  = cminmax>=para{1};
highbordermask = cminmax >= max(para{3},para{1});
if ~isempty(mask)
    bordermask(mask(:,:,1) > 0) = 1;
    bordermask(mask(:,:,2) > 0) = 0;
    bordermask(mask(:,:,3) > 0) = 0;
    highbordermask(mask(:,:,1) > 0) = 1;
    highbordermask(mask(:,:,2) > 0) = 0;
    highbordermask(mask(:,:,3) > 0) = 0;
end
if para{2} >1
    highbordermask = bwareaopen(highbordermask,para{2});
end
bordermask = imreconstruct(highbordermask, bordermask,8);
