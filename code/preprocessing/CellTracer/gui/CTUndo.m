function CTUndo()
h = gcf;
handles = guidata(h);
if strcmp(handles.undo.type,'input')
    handles.undo.flag = 0;
    handles.im{handles.io.index} = handles.undo.data;
    handles.undo.type = 'None';
    guidata(h,handles);
    showFrame(h,[]);
    ax1 = getAxes(gcf,1);subplot(ax1);
    temp = handles.im{handles.io.index}{1}.raw;
    imshow(temp);
    updateMenuStatus();
else
    msgbox('Undo not defined yet');
end
msgbox('To be finished');