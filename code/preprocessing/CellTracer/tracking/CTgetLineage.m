function result = CTgetLineage(h,frame)
handles = guidata(h);
proj = handles.project;
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end

filename = proj.filename;
outputpath = proj.outputpath;
currentpath = pwd;
cd(outputpath);
if iscell(filename)
    frames  = length(filename);
else
    frames = 1;
    filename = cellstr(filename);
end
matfile = [char(filename{frame}),'.lineage.mat'];
if exist(matfile,'file')
    load(matfile);
    result = lineage;
else
    result = [];
end
cd(currentpath);

