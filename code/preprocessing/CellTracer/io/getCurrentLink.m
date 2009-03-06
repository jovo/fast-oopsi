function [clink ccolor alignment] = getCurrentLink(h,frame)
clink = []; ccolor = [];
handles = guidata(h);
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end

%get what data to show
totalframes = handles.io.frames;

%populate link info if there is nothing there before
if isempty(handles.link) 
    for frame = 1: totalframes
        handles.link{frame}.corres = [];
        handles.link{frame}.colors = [];
        handles.link{frame}.alignment = [];
    end
end
if frame < totalframes && frame >= 1   
    clink = handles.link{frame}.corres;
    ccolor = handles.link{frame}.colors;
    try
        alignment = handles.link{frame}.alignment;
    catch
        alignment = [];
    end
else
    clink = [];
    ccolor = [];
    alignment = [];
end
guidata(h,handles);