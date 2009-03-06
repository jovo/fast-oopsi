function CTrestoreFromSnapshot(h,pathname,snapfilename)
handles = guidata(h);
currentpath = pwd;

project = handles.project;
filename = project.filename;
outputpath = project.outputpath;
if iscell(filename)
    frames  = length(filename);
else
    frames = 1;
    filename = cellstr(filename);
end
for frame = 1 : frames
    sourcefile = fullfile(outputpath, [char(filename{frame}),'.mat']);
    destfile = fullfile(pathname,'snapshot',[char(filename{frame}),'.mat']);
    copyfile(destfile,sourcefile);
end

cd(pathname);
load(snapfilename,'-mat');
cd(currentpath);

handles.viewtype = proj.viewtype;
handles.viewindex = proj.viewindex;
%handles.project = proj.project;
handles.io = proj.io;
handles.link = proj.link;
handles.undo = proj.undo; 
guidata(h,handles);