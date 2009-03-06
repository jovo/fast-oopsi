function updateCellButtons(v,frames,handles)
if v==1 
    set(handles.FirstCell,'Enable','off');
    set(handles.PreviousCell,'Enable','off');   
else
    set(handles.FirstCell,'Enable','on');
    set(handles.PreviousCell,'Enable','on');   
end
if v==frames
    set(handles.NextCell,'Enable','off');
    set(handles.LastCell,'Enable','off');
else
    set(handles.NextCell,'Enable','on');
    set(handles.LastCell,'Enable','on');
end