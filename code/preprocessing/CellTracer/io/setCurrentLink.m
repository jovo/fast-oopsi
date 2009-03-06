function  setCurrentLink(h,clink,ccolor,alignment,frame)
handles = guidata(h);
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end

%get what data to show
totalframes = handles.io.frames;
if frame < totalframes
    handles.link{frame}.corres = clink;
    handles.link{frame}.colors = ccolor;
    handles.link{frame}.alignment = alignment;
end
guidata(h,handles);
