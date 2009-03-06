function para = CTgetParameters()
h = gcf;
handles = guidata(h);

method = get(handles.ButtonRun,'TooltipString');
switch lower(method)
   case {'background1'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.backgroundspread = 0.1;
        else
            para.backgroundspread = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.disk = 0;
        else
            para.disk = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.padding = 1;
        else
            para.padding = str2num(v);
        end
   case {'background2'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.mincellvolume = 20;
        else
            para.mincellvolume = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.maxcellvolume = 20000;
        else
            para.maxcellvolume = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.disk = 0;
        else
            para.disk = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.minintensity = -1;
        else
            para.minintensity = str2num(v);
        end
        v = get(handles.editSecondary3,'String'); 
        if strcmp(v,'auto')
            para.maxintensity = -1;
        else
            para.maxintensity = str2num(v);
        end
        v = get(handles.editSecondary4,'String'); 
        if strcmp(v,'auto')
            para.maxseed = -1;
        else
            para.maxseed = str2num(v);
        end
    case {'background3'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.maxhalfwidth = 20;
        else
            para.maxhalfwidth = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.backgroundspread = 60;
        else
            para.backgroundspread = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.disk = 20000;
        else
            para.disk = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.greedy = 0;
        else
            para.greedy = str2num(v);
        end
        v = get(handles.editSecondary3,'String'); 
        if strcmp(v,'auto')
            para.fillholes = 0;
        else
            para.fillholes = str2num(v);
        end
    case {'border1'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.borderspread = 0.5;
        else
            para.borderspread = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.minsize = 1;
        else
            para.minsize = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.borderthreshold = 0.0;
        else
            para.borderthreshold = str2num(v);
        end
        v = get(handles.editSecondary3,'String'); 
        if strcmp(v,'auto')
            para.maxhalfwidth = 40;
        else
            para.maxhalfwidth = str2num(v);
        end
    case {'border2'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.borderspread = 0.5;
        else
            para.borderspread = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
         if strcmp(v,'auto')
            para.maxhalfwidth = 20;
        else
            para.maxhalfwidth = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.rankingmethod = 1;
        else
            para.rankingmethod = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.borderthreshold = 0.25;
        else
            para.borderthreshold = str2num(v);
        end
        v = get(handles.editSecondary3,'String'); 
        if strcmp(v,'auto')
            para.disk = 3;
        else
            para.disk = str2num(v);
        end
     case {'border3'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.borderspread = 0.5;
        else
            para.borderspread = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.minvolume = 20;
        else
            para.minvolume = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.rankingmethod = 1;
        else
            para.rankingmethod = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.borderthreshold = 0.25;
        else
            para.borderthreshold = str2num(v);
        end
        v = get(handles.editSecondary3,'String'); 
        if strcmp(v,'auto')
            para.maxhalfwidth = 20;
        else
            para.maxhalfwidth = str2num(v);
        end
        v = get(handles.editSecondary4,'String'); 
        if strcmp(v,'auto')
            para.disk = 3;
        else
            para.disk = str2num(v);
        end
    case {'border4'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.maxsize = 20;
        else
            para.maxsize = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.minpercentile = 0.75;
        else
            para.minpercentile = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.maxpercentile = 1;
        else
            para.maxpercentile = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.globalthreshold = 255;
        else
            para.globalthreshold = str2num(v);
        end
        v = get(handles.editSecondary3,'String'); 
        if strcmp(v,'auto')
            para.disk = 3;
        else
            para.disk = str2num(v);
        end
    case {'cell1'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.cellscore = 1.0;
        else
            para.cellscore = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.minvolume = 50;
        else
            para.minvolume = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
         if strcmp(v,'auto')
            para.disk = 0.25;
        else
            para.disk = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.dispercentile = 1;
        else
            para.dispercentile = str2num(v);
        end
    case {'cell2'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.maxsize = 1.0;
        else
            para.maxsize = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.percentile = 75;
        else
            para.percentile = str2num(v)*100;
        end
        v = get(handles.editPrimary3,'String'); 
         if strcmp(v,'auto')
            para.cellscore = 10.0;
        else
            para.cellscore = str2num(v);
         end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.disk = 1;
        else
            para.disk = str2num(v);
        end
        v = get(handles.editSecondary2,'String'); 
        if strcmp(v,'auto')
            para.celldeg = 1;
        else
            para.celldeg = str2num(v);
        end
    case {'globalalignment'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.shift = 20;
        else
            para.shift = str2num(v);
        end
    case {'celltracking'}
        v = get(handles.editPrimary1,'String'); 
        if strcmp(v,'auto')
            para.shift = 20;
        else
            para.shift = str2num(v);
        end
        v = get(handles.editPrimary2,'String'); 
        if strcmp(v,'auto')
            para.minscore = 0.1;
        else
            para.minscore = str2num(v);
        end
        v = get(handles.editPrimary3,'String'); 
        if strcmp(v,'auto')
            para.celldisplacement = 20;
        else
            para.celldisplacement = str2num(v);
        end
        v = get(handles.editSecondary1,'String'); 
        if strcmp(v,'auto')
            para.shiftscalefactor = 1;
        else
            para.shiftscalefactor = str2num(v);
        end
    otherwise
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
handles.parameters = para;
guidata(h,handles);

