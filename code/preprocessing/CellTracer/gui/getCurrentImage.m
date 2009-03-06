function [im mask] = getCurrentImage(h,frame)
im = []; mask = [];
handles = guidata(h);
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end
viewtype = handles.viewtype; %get which type of data to show 
if strcmp(viewtype,'input') 
    %get what data to show
    viewindex = handles.viewindex;
    totalframes = handles.io.frames;
    if frame > totalframes || frame < 1  %out of boundary
        [im mask] = CTreadImage(handles.project,1,viewindex);
        im(im>0) = 0;
        mask= [];
    else
        [im mask] = CTreadImage(handles.project,frame,viewindex);
    end
    if isempty(mask)
        mask(:,:,1) = im<0;
        mask(:,:,2) = mask(:,:,1);
        mask(:,:,3) = mask(:,:,1);
    end
elseif  strcmp(viewtype,'segmentation') %segmentation
    %get what data to show
    viewindex = handles.viewindex;
    totalframes = handles.io.frames;
    if frame > totalframes || frame < 1  %out of boundary
        [im mask] = CTreadImage(handles.project,1,viewindex);
        im(im>0) = 0;
        mask= [];
    else
        [im mask] = CTreadImage(handles.project,frame,viewindex);
    end
    if isempty(mask)
        mask(:,:,1) = im<0;
        mask(:,:,2) = mask(:,:,1);
        mask(:,:,3) = mask(:,:,1);
    else
        L = bwlabel(mask(:,:,3) > 0);
        mask = label2rgb(L, 'spring', 'c', 'shuffle'); 
    end
    
elseif  strcmp(viewtype,'tracking') %tracking
    temp = handles.viewindex;
    handles.viewindex = -1;
    handles.viewtype = 'input'; 
    guidata(h,handles);
    [dummy mask1] = getCurrentImage(h,frame);
    [dummy mask2] = getCurrentImage(h,frame+1);
    im = uint8(mask1(:,:,3) > 0) * 255;         %left panel is always 256 grey-scale
    mask = mask2(:,:,3) > 0;
    handles.viewtype = viewtype;
    handles.viewindex = temp;
    guidata(h,handles);
else
     msgbox('Not defined');
end