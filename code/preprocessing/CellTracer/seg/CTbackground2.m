function backgroundmask = CTbackground2(im,mask,para,debug)
    im1 = im; 
    minvolume = para{1};
    maxvolume = para{2};
    disk = para{3};
    minintensity = para{4};
    maxintensity = para{5};
    maxseed = para{6};
    limit = [minintensity minintensity maxseed maxintensity];
    if ~isempty(mask)
        im1(mask(:,:,1) > 0) = 255;
        im1(mask(:,:,2) > 0) = 255;
        [a b] = CTScreen(im1,mask(:,:,3),maxvolume,minvolume,disk,8,limit,[],[],[]);
    else
        [a b] = CTScreen(im1,[],maxvolume,minvolume,disk,8,limit,[],[],[]);
    end
    
if debug>0
    figure,imshow(a);title('result 0');
    figure,imshow(b);title('result 1');
end

datamask = a > 0;

datamask = CTiterativeApply(datamask,0,4,{@imfill},'holes');
backgroundmask = datamask ==0; %to get the background mask
if ~isempty(mask)
    backgroundmask(mask(:,:,2) > 0) = 1;
    backgroundmask(mask(:,:,1) > 0) = 0;
    backgroundmask(mask(:,:,3) > 0) = 0;
end
