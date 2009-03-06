function result = showFrame(hObject,frame)
    cache = [];
    result = 1;
    handles = guidata(hObject);
    viewtype = handles.viewtype;
    if ~isempty(handles.project)
        set(handles.StatusBar,'String', 'Switching between frames.... Please Wait...'); pause(0.1);
        if ~isempty(frame) && frame == -1
            frame = [];
            rescale = 0;
        else
            rescale  = 1;
        end
        try
            [im mask] = getCurrentImage(hObject,frame);    
            if ~strcmp(viewtype,'tracking')     
                ax2 = getAxes(gcf,2); v = axis; subplot(ax2); 
                if max(mask(:)) <=1
                    imshow(uint8(mask) * 255);
                else
                    imshow(mask);
                end
                ax1 = getAxes(gcf,1);subplot(ax1); imshow(im); 
                if rescale  >0 
                    axis(v);
                else
                    axis image;
                    linkaxes([ax1 ax2]);
                end
            else
                [clink ccolor] = getCurrentLink(hObject,[]);  %get all verified links
                if ~isempty(clink) 
                    ccolor = CTassignLinkColor(clink, ccolor);
                    set(handles.StatusBar,'String', 'Regenerating lineage images.... Please Wait...'); pause(0.1);
                    [left right] = redrawImages(im,mask,ccolor);
                    cache.left = left;
                    cache.right = right;
                    cache.frame = frame;
                    handles.cache = cache;
                    guidata(hObject,handles);
                else
                    ax2 = getAxes(gcf,2); v = axis; subplot(ax2); 
                    if max(mask(:)) <=1
                        imshow(uint8(mask) * 255);
                    else
                        imshow(mask);
                    end
                    ax1 = getAxes(gcf,1);subplot(ax1); imshow(im); 
                    if rescale  >0 
                        axis(v);
                    else
                        axis image;
                        linkaxes([ax1 ax2]);
                    end
                end
            end
        catch
            result = 0;
        end
        set(handles.StatusBar,'String', 'Switch between images.... Done');
    end
end
