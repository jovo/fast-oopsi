function CTsetParameters(method)
h = gcf;
handles = guidata(h);
set(handles.ButtonRun,'TooltipString',method);
handles.parameters.method = method;
if strcmp(handles.viewtype,'input')
    set(handles.checkboxAllFrames,'Enable','On');
    set(handles.checkboxDebug,'Enable','On');
    set(handles.checkboxSequential,'Enable','On');
    set(handles.ButtonRun,'Enable','On');
elseif strcmp(handles.viewtype,'tracking')
    set(handles.checkboxAllFrames,'Enable','On');
    set(handles.checkboxDebug,'Enable','On');
    set(handles.checkboxSequential,'Enable','Off');
    set(handles.ButtonRun,'Enable','On');
else
    set(handles.checkboxAllFrames,'Enable','Off');
    set(handles.checkboxDebug,'Enable','Off');
    set(handles.checkboxSequential,'Enable','Off');
    set(handles.ButtonRun,'Enable','Off');
end

if (get(handles.checkboxAllFrames,'Value') == get(handles.checkboxAllFrames,'Max'))
    para.allframes = 1;
else
	para.allframes = 0;
end
if (get(handles.checkboxDebug,'Value') == get(handles.checkboxDebug,'Max'))
    para.debug = 1;
else
	para.debug = 0;
end

if (get(handles.checkboxSequential,'Value') == get(handles.checkboxSequential,'Max'))
    para.sequential = 1;
else
	para.sequential = 0;
end

allcontrols = get(handles.uipanelPrimary,'Children');
for i = 1:length(allcontrols)
    set(allcontrols(i),'Visible','Off');
end
allcontrols = get(handles.uipanelSecondary,'Children');
for i = 1:length(allcontrols)
    set(allcontrols(i),'Visible','Off');
end

switch lower(method)
   case {'background1'}
       enableParameter(handles,1,1,'Background Spread:','0.1');
       enableParameter(handles,2,1,'Structure Element Radius:','0');
       enableParameter(handles,2,2,'Padding Method (0 or 1):','1');
   case {'background2'}
       enableParameter(handles,1,1,'Minimum Cell Volume:','20');
       enableParameter(handles,1,2,'Maximum Cell Volume:','20000');
       
       enableParameter(handles,2,1,'Structure Element Radius:','auto');
       enableParameter(handles,2,2,'Minimum Cell Intensity:','auto');
       enableParameter(handles,2,3,'Maximum Cell Intensity:','auto');
       enableParameter(handles,2,4,'Maximum Cell Intensity Seed:','auto');
   case {'background3'}
       enableParameter(handles,1,1,'Maximum Cell Half Width:','20');
       enableParameter(handles,1,2,'Background Intensity spread:','60');
       enableParameter(handles,2,1,'Structure Element Radius:','0');
       enableParameter(handles,2,2,'Greedy:','0');
       enableParameter(handles,2,3,'Fill Holes:','0');
   case {'border1'}
       enableParameter(handles,1,1,'Minimum Ranking Threshold:','0.5');
       enableParameter(handles,2,1,'Minimum Border Volume:','auto');
       enableParameter(handles,2,2,'Upper Ranking Threshold:','auto');
       enableParameter(handles,2,3,'Maximum Cell Half Width:','auto');
   case {'border2'}
       enableParameter(handles,1,1,'Minimum Ranking Threshold:','0.5');
       enableParameter(handles,1,2,'Maximum Cell Half Width:','20');
       
       enableParameter(handles,2,1,'Ranking Method(1 or 3):','1');
       enableParameter(handles,2,2,'Count Transform Threshold:','auto');
       enableParameter(handles,2,3,'Structure Element Radius:','auto');
   case {'border3'}
       enableParameter(handles,1,1,'Minimum Ranking Threshold:','0.5');
       enableParameter(handles,1,2,'Minimum Border Volume:','20');
       
       enableParameter(handles,2,1,'Ranking Method(1 or 3):','1');
       enableParameter(handles,2,2,'Count Transform Threshold:','0.25');
       enableParameter(handles,2,3,'Maximum Cell Half Width:','20');
       enableParameter(handles,2,4,'Minimum Cell Half Width:','3');
   case {'border4'}
       enableParameter(handles,1,1,'Maximum Cell Half Width:','20');
       enableParameter(handles,1,2,'Lower Ranking Threshold:','0.75');
       
       enableParameter(handles,2,1,'Upper Ranking Threshold:','1.0');
       enableParameter(handles,2,2,'Global Intensity Threshold:','255');
       enableParameter(handles,2,3,'Structure Element Radius:','3');
   case {'cell1'}
       enableParameter(handles,1,1,'Maximum Cell Score:','1.0');
       enableParameter(handles,1,2,'Minimum Cell Volume:','50');
       
       enableParameter(handles,2,1,'Structure Element Radius:','1');
       enableParameter(handles,2,2,'Smoothing Parameter:','0');
   case {'cell2'}
       enableParameter(handles,1,1,'Maximum Cell Half Width:','20');
       enableParameter(handles,1,2,'Intensity threshold:','0.75');
       enableParameter(handles,1,3,'Cutoff score:','10');
       
       enableParameter(handles,2,1,'Structure Element Radius:','1');
       enableParameter(handles,2,2,'Curve degree :','1');
    case {'globalalignment'}
       enableParameter(handles,1,1,'Maximum Shift:','20');
    case {'celltracking'}
       enableParameter(handles,1,1,'Neighborhood Size:','20');
       enableParameter(handles,1,2,'Minimum Overlapping Score:','0.1');
       enableParameter(handles,1,3,'Maxinum Cell Displacement:','20');
       
       enableParameter(handles,2,1,'Neighborhood Scale Factor:','1');
   otherwise
end

function enableParameter(handles,type,which, str,defaultvalue)
if type == 1 %primary
    allcontrols = get(handles.uipanelPrimary,'Children');
    switch which
       case {1}
           ch = findall(allcontrols,'Tag','textPrimary1');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editPrimary1');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
       case {2}
           ch = findall(allcontrols,'Tag','textPrimary2');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editPrimary2');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
       case {3}
           ch = findall(allcontrols,'Tag','textPrimary3');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editPrimary3');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
        otherwise
            msgbox('Too many parameters');
    end
else
    allcontrols = get(handles.uipanelSecondary,'Children');
    switch which
       case {1}
           ch = findall(allcontrols,'Tag','textSecondary1');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editSecondary1');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
       case {2}
           ch = findall(allcontrols,'Tag','textSecondary2');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editSecondary2');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
       case {3}
           ch = findall(allcontrols,'Tag','textSecondary3');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editSecondary3');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
       case {4}
           ch = findall(allcontrols,'Tag','textSecondary4');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editSecondary4');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
       case {5}
           ch = findall(allcontrols,'Tag','textSecondary5');
           set(ch,'Visible','On');
           set(ch,'String',str);
           ch = findall(allcontrols,'Tag','editSecondary5');
           set(ch,'Visible','On');
           set(ch,'Enable','On');
           set(ch,'String',defaultvalue);
        otherwise
            msgbox('Too many parameters');
    end
end