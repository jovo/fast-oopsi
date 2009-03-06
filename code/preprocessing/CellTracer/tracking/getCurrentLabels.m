function result = getCurrentLabels(h,frame)
handles = guidata(h);
result = [];
if ~strcmp(handles.viewtype,'tracking')
    msgbox('getCurrentLabels::program bug');
else
    if isempty(frame) || frame < 1
        frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
    end
    try
        %if there is already labels, then return it
        result = handles.link{frame}.labels;
    catch
        %no labels yet, so create
        [im mask] = getCurrentImage(h,frame);
        [result num] = bwlabel(im>0,4); 
        if num <=255
            handles.link{frame}.labels = uint8(result);
        else
            handles.link{frame}.labels = uint16(result);
        end
        guidata(h,handles);
    end
    
end
