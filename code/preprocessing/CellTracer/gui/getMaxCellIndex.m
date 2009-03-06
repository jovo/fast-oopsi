function result = getMaxCellIndex(h, frame)
handles = guidata(h);
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end
bwlabels = getCurrentLabels(h,frame);
result = max(bwlabels(:));