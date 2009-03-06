function  purgeCurrentImage(h,frame,index)
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
matfile = [char(filename{frame}),'.mat'];
load(matfile);
totalviews = length(views);
if totalviews > index
    for i =1:index
        newviews{i} = views{i};
    end
    views = newviews;
    save(matfile, 'views');
end
cd(currentpath);
guidata(h,handles);
