function showCell(hObject,cell)
    handles = guidata(hObject);
    %viewtype = handles.viewtype; %viewtype should always be 'tracking'
    if ~isempty(handles.project)
        %deal with the left panel first
        frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
        [im mask] = getCurrentImage(hObject,frame);
        bwlabels = getCurrentLabels(hObject,frame);
        [clink ccolor] = getCurrentLink(hObject,frame);  %get all verified links
        if ~isempty(clink) 
            try
                if handles.cache.frame ~= frame
                    msgbox('Cached lineage images need to be updated!');
                    handles.cache=[];
                    temp = handles.cache.left; %trigger a catch
                end
                leftimage = handles.cache.left;
                bw = bwlabels == cell; 
                if any(bw(:))
                    bw1 = bwperim(bw);
                    bw(bw1>0) = 0;
                    bw2 = bwperim(bw);
                   
                end
                for i = 1: 3
                    temp = leftimage(:,:,i);
                    temp(bw1 > 0 | bw2 > 0) = 255;
                    leftimage(:,:,i) = temp;
                end
                ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
                imshow(leftimage);
            catch
                ccolor = CTassignLinkColor(clink, ccolor);
                bw = bwlabels == cell;
                if any(bw(:))
                    bw1 = bwperim(bw);
                    bwlabels(bw1>0) = 0;
                    bw(bw1>0) = 0;
                    bw1 = bwperim(bw);
                    bwlabels(bw1>0) = 0;
                end
                %redraw the left axis based on view
                ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
                imshow(drawimage(3,1,im,ccolor,bwlabels));
            end
            
            try 
                rightimage = handles.cache.right;
                bwmask =  getCurrentLabels(hObject,frame+1);
                if isempty(bwmask) 
                    bwmask = mask;
                    bwmask = bwlabel(bwmask,4); 
                    %redraw the right axis based on view
                    ax2 = getAxes(gcf,2); subplot(ax2);
                    imshow(drawimage(3,2,mask,ccolor,bwmask));
                else
                    linklist = find(clink(cell,:) >0);
                    if ~isempty(linklist)
                        for i = 1:length(linklist)
                            bw = bwmask == linklist(i);
                            if any(bw(:))
                                if any(bw(:))
                                    bw1 = bwperim(bw);
                                    bw(bw1>0) = 0;
                                    bw2 = bwperim(bw);
                                end
                                for i = 1: 3
                                    temp = rightimage(:,:,i);
                                    temp(bw1 > 0 | bw2 > 0) = 255;
                                    rightimage(:,:,i) = temp;
                                end
                            end
                        end
                    end
                    %redraw the right axis based on view
                    ax2 = getAxes(gcf,2); subplot(ax2);
                    imshow(rightimage);
                end
            catch
                bwmask =  getCurrentLabels(hObject,frame+1);
                if isempty(bwmask) 
                    bwmask = mask;
                    bwmask = bwlabel(bwmask,4); 
                else
                    linklist = find(clink(cell,:) >0);
                    if ~isempty(linklist)
                        for i = 1:length(linklist)
                            bw = bwmask == linklist(i);
                            if any(bw(:))
                                bw1 = bwperim(bw);
                                bwmask(bw1>0) = 0;
                                bw(bw1>0) = 0;
                                bw1 = bwperim(bw);
                                bwmask(bw1>0) = 0;
                            end
                        end
                    end
                end
                %redraw the right axis based on view
                ax2 = getAxes(gcf,2); subplot(ax2);
                imshow(drawimage(3,2,mask,ccolor,bwmask));
            end
            
            
            axis(v);
        else
            ax2 = getAxes(gcf,2); v = axis; subplot(ax2); 
            if max(mask(:)) <=1
                imshow(uint8(mask) * 255);
            else
                imshow(mask);
            end
            
            bw = bwlabels == cell;
            if any(bw(:))
                im(bw>0) = 0;
                bw1 = bwperim(bw);
                im(bw1>0) = 1; 
                bw(bw1>0) = 0;
                bw1 = bwperim(bw);
                im(bw1>0) = 1;
                im(im>0) = 255;
            end
            ax1 = getAxes(gcf,1);subplot(ax1); imshow(im); 
            axis(v);
        end
    end
end
