function CTcacheLineageImages(hObject,frame)
handles = guidata(hObject);
if strcmp(handles.viewtype,'tracking')
    if ~isempty(handles.project)
       try
           set(handles.StatusBar,'String', 'Caching lineage images.... Plelase wait');
           pause(0.1);
           if isempty(frame)
                frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
           end
            [clink ccolor] = getCurrentLink(hObject,frame);  %get all verified links
            if ~isempty(clink)
                bwlabels = getCurrentLabels(hObject,frame);
                im = bwlabels > 0;
            
             
                ccolor = CTassignLinkColor(clink, ccolor);
                handles.cache.left = drawimage(3,1,im,ccolor,bwlabels);
            
                bwlabels = getCurrentLabels(hObject,frame+1);
                im = bwlabels > 0;
            
                ccolor = CTassignLinkColor(clink, ccolor);
                handles.cache.right = drawimage(3,2,im,ccolor,bwlabels);
                handles.cache.frame = frame;
                guidata(hObject,handles);
            end
            set(handles.StatusBar,'String', 'Caching lineage images... Done');
       catch 
           set(handles.StatusBar,'String', 'Caching lineage images... Failed');
       end
    end
else
    msgbox('Please switch to Tracking View beofre running this tool'); 
end
