function setCurrentLabels(h,labels,frame)
handles = guidata(h);
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end
handles.link{frame}.labels = uint8(labels);
guidata(h,handles);
