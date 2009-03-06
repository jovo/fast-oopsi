function updateButtons(v,frames,handles)
if v==1 
    set(handles.FirstButton,'Enable','off');
    set(handles.PreviousButton,'Enable','off');   
else
    set(handles.FirstButton,'Enable','on');
    set(handles.PreviousButton,'Enable','on');   
end
if v==frames
    set(handles.NextButton,'Enable','off');
    set(handles.LastButton,'Enable','off');
else
    set(handles.NextButton,'Enable','on');
    set(handles.LastButton,'Enable','on');
end