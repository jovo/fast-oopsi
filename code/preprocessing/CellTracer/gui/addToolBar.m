function addToolBar()
    fh = findobj('Tag','figure1');
    set(fh,'Toolbar','figure');    % unhide Hide the standard toolbar
    tbh = findall(fh,'Type','uitoolbar');
    oldorder = allchild(tbh);
    %delete(oldorder([1 2 3 4 6 10 11 12 13 14]));
    delete(oldorder([1 2 3 4 6 10 12 13 14]));
%     dch = findall(tbh,'TooltipString','Data Cursor');
%     if ~isempty(dch)
%         get(dch,'Type')
%     end
    load icons.mat icons;
    % Add a toggle tool to the toolbar
    %-------------------------------------------------------  
    pth6 = uitoggletool(tbh,'CData',icons{6},...
           'TooltipString','Overlay Tool',...
           'HandleVisibility','on',...
           'onCallback',@oncrop_callback,...
           'offcallback',@offcrop_callback);
       function oncrop_callback(hObject,eventdata)
           updatemask(hObject,0);
       end  
       function offcrop_callback(hObject,eventdata)
           updatemask(hObject,-1);
       end
    
    pth4 = uipushtool(tbh,'CData',icons{4},...
           'TooltipString','Border Region Mask Tool',...
           'HandleVisibility','on',...
           'Tag','BorderMask',...
           'ClickedCallback',@border_callback);
       function border_callback(hObject,eventdata)
           updatemask(hObject,1);
       end  
    
    pth2 = uipushtool(tbh,'CData',icons{2},...
               'TooltipString','Background Region Mask Tool',...
               'HandleVisibility','on',...
               'Tag','BackgroundMask',...
               'ClickedCallback',@backgroundmask_callback);
           function backgroundmask_callback(hObject,eventdata)
               updatemask(hObject,2);
           end
    pth3 = uipushtool(tbh,'CData',icons{3},...
           'TooltipString','Cell Object Mask Tool',...
           'HandleVisibility','on',...
           'Tag','CellMask',...
           'ClickedCallback',@cell_callback);
       function cell_callback(hObject,eventdata)
           updatemask(hObject,3);
       end
    pth5 = uipushtool(tbh,'CData',icons{5},...
           'TooltipString','Object Select Tool',...
           'HandleVisibility','on',...
           'ClickedCallback',@object_callback);
       function object_callback(hObject,eventdata)
           selectobject(hObject);
       end
    function selectobject(hObject)
        handles = guidata(hObject);
        frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
        
        [im mask] = getCurrentImage(hObject,[]);  %get curent images based on view
        [clink ccolor] = getCurrentLink(hObject,[]);  %get all verified links
        currentview = getCurrentView(hObject);
        
        %label the left axis
        bwim = im >0;
        if currentview == 3
            bwim = getCurrentLabels(hObject,frame);
            numim = max(bwim(:));
        else
            [bwim numim] = bwlabel(bwim,4); 
        end
        objectselected{1} = zeros(numim,3);
        
        %label the right panel
        if length(size(mask)) == 3
            [bwmask1 nummask1] = bwlabel(mask(:,:,1),4); 
            [bwmask2 nummask2] = bwlabel(mask(:,:,2),4);
            [bwmask3 nummask3] = bwlabel(mask(:,:,3),4);
            nummask = nummask1 + nummask2 + nummask3;
            bwmask(:,:,1) = bwmask1;
            bwmask(:,:,2) = bwmask2;
            bwmask(:,:,3) = bwmask3;
        else
            if currentview == 3
                bwmask = getCurrentLabels(hObject,frame+1);
                if isempty(bwmask)
                    bwmask = mask;
                    [bwmask nummask] = bwlabel(bwmask,4); 
                else
                    nummask = max(bwmask(:));
                end
            else
                bwmask = mask;
                [bwmask nummask] = bwlabel(bwmask,4); 
            end
        end
        objectselected{2} = zeros(nummask,3);
        
        %redraw the left axis based on view
        ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
        imshow(drawimage(currentview,1,im,ccolor,bwim));
        
        %redraw the right axis based on view
        ax2 = getAxes(gcf,2); subplot(ax2);
        imshow(drawimage(currentview,2,mask,ccolor,bwmask));
        axis(v);
        
        [xi,yi,but] = ginput(1);
        while but == 1
            xi = uint16(xi); yi = uint16(yi);
            if xi ==0, xi = 1; end
            if yi ==0, yi = 1; end
            if gca == ax1
                objectselected = updateSelectedObjects(currentview, 1, objectselected,xi, yi,ccolor,v,bwim,im);
            elseif gca == ax2
                objectselected = updateSelectedObjects(currentview, 2, objectselected,xi, yi,ccolor,v,bwmask,mask);
            end
            [xi,yi,but] = ginput(1);
        end     %end of selectobject
        
        
        function updatedObjects = updateSelectedObjects(currentview, whichaxis, objectselected,xi, yi,ccolor,v,labels,img)
            maxobjects{1} = 1; maxobjects{2} = 2; %used only in tracking view to restrict number of objects selected
            try 
                if length(size(labels)) == 3
                    n1 = max(labels(:,:,1));
                    n2 = max(labels(:,:,2));
                    if labels(yi,xi,1) > 0
                        currentobject = labels(yi,xi,1);
                    elseif labels(yi,xi,2) > 0
                        currentobject = labels(yi,xi,2) + n1;
                    elseif labels(yi,xi,3) > 0
                        currentobject = labels(yi,xi,3) + n1 + n2;
                    else
                        currentobject = 0;
                    end
                else
                    currentobject = labels(yi,xi);
                end
                if currentobject >0
                    if currentview <3 %non tracking view, allow choosing many objects
                        if objectselected{whichaxis}(currentobject,1) > 0
                            objectselected{whichaxis}(currentobject,1:3) = 0;
                        else
                            objectselected{whichaxis}(currentobject,1) = 1;
                            objectselected{whichaxis}(currentobject,2) = xi;
                            objectselected{whichaxis}(currentobject,3) = yi;
                        end
                    else
                        [pm pn] = size(objectselected{whichaxis});
                        if objectselected{whichaxis}(currentobject,1) > 0
                            temp = objectselected{whichaxis}(currentobject,1);
                            for i = 1: pm
                                if objectselected{whichaxis}(i,1) > temp
                                    objectselected{whichaxis}(i,1) = objectselected{whichaxis}(i,1) - 1;
                                end
                            end
                            objectselected{whichaxis}(currentobject,1:3) = 0;
                        else
                            numselected = length(find(objectselected{whichaxis}(:,1)>0));
                            if numselected < maxobjects{whichaxis}
                                objectselected{whichaxis}(currentobject,1) = numselected + 1;
                                objectselected{whichaxis}(currentobject,2) = xi;
                                objectselected{whichaxis}(currentobject,3) = yi;
                            else
                                for i = 1: pm
                                    if objectselected{whichaxis}(i,1) > 0
                                        objectselected{whichaxis}(i,1) = objectselected{whichaxis}(i,1) - 1;
                                        if objectselected{whichaxis}(i,1) == 0
                                            objectselected{whichaxis}(i,1:3) = 0;
                                        end
                                    end
                                end 
                                objectselected{whichaxis}(currentobject,1) = numselected;
                                objectselected{whichaxis}(currentobject,2) = xi;
                                objectselected{whichaxis}(currentobject,3) = yi;
                            end
                        end
                    end
                    
                    imshow(drawimage(currentview,whichaxis,img,ccolor,labels)); axis(v);
                    [pm pn] = size(objectselected{whichaxis});
                    hold on
                    for i = 1:pm
                        if objectselected{whichaxis}(i,1) > 0;
                            plot(objectselected{whichaxis}(i,2),objectselected{whichaxis}(i,3),'c+','UIContextMenu',addcmenu(whichaxis,objectselected,currentview));
                        end
                    end
                    hold off
                end
            catch
                
            end
            updatedObjects = objectselected;
        end
        function cmh = addcmenu(whichaxes,points,currentview)
            cmh = uicontextmenu;
            if whichaxes == 1
                if currentview ==1
                   %do nothing for now
                elseif currentview ==2
                    %do nothing for now
                elseif currentview ==3
                    cb1 = 'CellTracer(''cmenu_linkmaskedobject'',gcbo,[],guidata(gcbo))';
                    cb2 = 'CellTracer(''cmenu_removemaskedobject'',gcbo,[],guidata(gcbo))';
                    item1 = uimenu(cmh, 'Label', 'Link/Unlink selected blob(s)', 'Callback', cb1,'UserData',points);
                    item2 = uimenu(cmh, 'Label', 'Remove selected blob(s)', 'Callback', cb2,'UserData',{points,currentview});
                end
            else
                if currentview ==1
                    cb1 = 'CellTracer(''cmenu_removemaskedobject'',gcbo,[],guidata(gcbo))';
                    item1 = uimenu(cmh, 'Label', 'Remove selected blob(s)', 'Callback', cb1,'UserData',{points,currentview});
                elseif currentview ==2
                    cb1 = 'CellTracer(''cmenu_removemaskedobject'',gcbo,[],guidata(gcbo))';
                    item1 = uimenu(cmh, 'Label', 'Remove selected blob(s)', 'Callback', cb1,'UserData',{points,currentview});
                elseif currentview ==3
                    cb1 = 'CellTracer(''cmenu_linkmaskedobject'',gcbo,[],guidata(gcbo))';
                    cb2 = 'CellTracer(''cmenu_removemaskedobject'',gcbo,[],guidata(gcbo))';
                    item1 = uimenu(cmh, 'Label', 'Link/Unlink selected blob(s)', 'Callback', cb1,'UserData',points);
                    item2 = uimenu(cmh, 'Label', 'Remove selected blob(s)', 'Callback', cb2,'UserData',{points,currentview});
                end
            end
        end        
    end
    function updatemask(hObject,which)   
        state = uisuspend(gcf);
        [im mask] = getCurrentImage(hObject,[]);
        handles = guidata(hObject);
        if which == -1
            ax1 = getAxes(gcf,1);subplot(ax1);v = axis; 
            imshow(im); axis(v);
        elseif which == 0
            ax1 = getAxes(gcf,1);subplot(ax1);v = axis;
            if ~strcmp(handles.viewtype,'tracking')
                for i = 1:3
                    temp = im; 
                    temp(mask(:,:,i)>0) = 128;
                    rgb(:,:,i) = temp;
                end 
            else
                mask = CTshiftImage(hObject,mask,[]);
                rgb(:,:,1) = uint8(im) * 255;
                rgb(:,:,2) = uint8(mask) * 255;
                rgb(:,:,3) = 0;
            end
            imshow(rgb); axis(v);
        else
            if ~strcmp(handles.viewtype,'tracking')
                ax1 = getAxes(gcf,1);subplot(ax1);v = axis;
                for i = 1:3
                    temp = im; temp(mask(:,:,i)>0) = 128;
                    rgb(:,:,i) = temp;
                end 
                imshow(rgb); axis(v);
                try 
                    bw = roipoly;
                catch
                    bw = im < 0;
                    uirestore(state);
                end
                if strcmp(handles.viewtype,'segmentation')
                    handles.viewtype = 'input';
                    guidata(hObject,handles);
                    [im1 mask] = getCurrentImage(hObject,[]);
                    handles.viewtype = 'segmentation';
                    guidata(hObject,handles);
                end
                for i = 1:3
                    temp = mask(:,:,i);
                    if i==which
                        temp = temp>0 | bw >0;
                    else
                        temp = temp > 0 & bw == 0;
                    end
                    mask(:,:,i) = temp;
                end
                try 
                    setCurrentImage(hObject,im,mask,[]);
                    ax2 = getAxes(gcf,2); subplot(ax2);
                    imshow(uint8(mask) * 255);
               
                    for i = 1:3
                        temp = im; temp(mask(:,:,i)>0) = 128;
                        rgb(:,:,i) = temp;
                    end 
                    subplot(ax1);
                    imshow(rgb);
                    axis(v);
                catch
                    uirestore(state);
                end

            else
                msgbox('Not defined yet');
            end
        end
        uirestore(state);
    end    
%     function crop(hObject)        
%         [im mask] = getCurrentImage(hObject,[]);
%         handles = guidata(hObject);
%         if strcmp(handles.viewtype,'input')
%             ax1 = getAxes(1);subplot(ax1);v = axis;
%             [A rect] = imcrop;
%             rect = uint16(rect); rect(rect<1) = 1;
%             
%             axis(v);
%         else
%             msgbox('Not defined yet');
%         end
%     end    
end
