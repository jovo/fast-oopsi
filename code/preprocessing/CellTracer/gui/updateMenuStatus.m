function updateMenuStatus()
h = gcf;
handles = guidata(h);

%update edit menu
hm = findall(h,'Tag','menuUndo');
%try 
    if handles.undo.flag > 0
        set(hm,'Enable','On');
    else
        set(hm,'Enable','Off');
    end
%end

%update view menu
ht = findall(h,'Tag','ViewTracking');
hs = findall(h,'Tag','ViewSegmentation');
hi = findall(h,'Tag','ViewInput');

%update toolbar
hborder = findall(h,'Tag','BorderMask');
hbackground = findall(h,'Tag','BackgroundMask');
hcell = findall(h,'Tag','CellMask');

%a = findall(tbh,'TooltipString','Data Cursor')
if strcmp(handles.viewtype,'input')
    set(ht,'Checked','Off');
    set(hs,'Checked','Off');
    set(hi,'Checked','On');
    set(hborder,'Enable','On');
    set(hbackground,'Enable','On');
    set(hcell,'Enable','On');
elseif strcmp(handles.viewtype,'segmentation')
    set(ht,'Checked','Off');
    set(hs,'Checked','On');
    set(hi,'Checked','Off');
    set(hborder,'Enable','On');
    set(hbackground,'Enable','On');
    set(hcell,'Enable','On');
else
    set(ht,'Checked','On');
    set(hs,'Checked','Off');
    set(hi,'Checked','Off');
    set(hborder,'Enable','Off');
    set(hbackground,'Enable','Off');
    set(hcell,'Enable','Off');
end

%cell control panel
allcontrols = get(handles.uipanelCell,'Children');
if strcmp(handles.viewtype,'tracking')
    for i = 1:length(allcontrols)
        set(allcontrols(i),'Visible','On');
    end
    set(handles.uipanelCell,'Visible','On');
else
    for i = 1:length(allcontrols)
        set(allcontrols(i),'Visible','Off');
    end
    set(handles.uipanelCell,'Visible','Off');
end