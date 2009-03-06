function result = CTiterativeSelectiveSegmentation2(im,inputmask,para,debug)
    maxsize = para{1};
    minpercentile = para{2};
    disk = para{3};
    cellscore = para{4};
    celldeg = para{5};

    if ~isempty(inputmask)
        im(inputmask(:,:,1) > 0) = 0;
        im(inputmask(:,:,2) > 0) = 0;
        im(inputmask(:,:,3) > 0) = 0;
    end
    mask = im > 0; L = bwlabel(mask,4); s = regionprops(L, 'Image','BoundingBox');  
    for k = 1:numel(s)
        bb = uint16(s(k).BoundingBox);
        a = im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1); img = a; 
        a(~s(k).Image) = 0; img(a>0) = 0;
        flag = getscore(cellscore,celldeg,a);
        if flag > cellscore
            b = a;

            estat = CTblock(b,[],[],maxsize, {@mean},[]);
            Mean = double(b) - estat(:,:,1);
            Mean = Mean - min(Mean(:));
            r1 = CTblock(Mean,[],[],maxsize,{@CTranks},{1});
            temp = r1 < minpercentile/100;
            temp = bwareaopen(temp,maxsize);
            if disk >0
                smoothed = CTiterativeOpen(temp, disk,1);
                smoothed = CTdistanceSmooth(smoothed,disk,0,debug);
                if disk > 2
                    smoothed = CTiterativeOpen(smoothed, 1,0);
                end
            else
                smoothed = temp;
            end

            b(a>0 & smoothed == 0) = 0;
            %restore
            c = b>0;  c = uint8(bwmorph(c,'thicken',disk * 2)); 
            a(a>0 & c==0) = 0;  
            img = img+a; 
            im(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1) = img;
        end
    end
    
    figure(101),imshow(im);
    result = uint8(im > 0);
    if ~isempty(inputmask)
        result(inputmask(:,:,3) >0) = 1;
    end
end


function score = getscore(cellscore,celldeg,a)
    if cellscore > 0
        tempbw = a==0;
        d =  bwdist(tempbw); temp = d(d>0);
        meand = mean(temp); std = sqrt(var(temp));
        d(d<=meand+std) = 0;  d(d>0) = 1;
        if length(find(d>0)) <= (celldeg+1)
            score = cellscore  + 1;
        else
            score = GetFlag(d,celldeg);
        end
    else
        score = cellscore + 1;
    end
end