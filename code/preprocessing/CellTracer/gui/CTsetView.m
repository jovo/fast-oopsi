function CTsetView(viewtype,viewindex)
handles = guidata(gcf);
handles.viewtype = viewtype;
handles.viewindex = viewindex;
guidata(gcf,handles);
