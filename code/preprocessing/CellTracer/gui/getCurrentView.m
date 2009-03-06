function cview = getCurrentView(hObject)
    handles = guidata(hObject); 
    if strcmp(handles.viewtype,'tracking')
        cview = 3;
    elseif strcmp(handles.viewtype,'segmentation')
        cview = 2;
    else %input
        cview = 1;
    end
end