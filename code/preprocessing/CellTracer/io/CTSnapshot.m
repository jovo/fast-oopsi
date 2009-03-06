function CTSnapshot(h,pathname,snapfilename)
handles = guidata(h);
currentpath = pwd;
if ~exist(pathname,'dir')
    mkdir(pathname);
end
cd(pathname);
try
    mkdir('snapshot');
catch
end

project = handles.project;
filename = project.filename;
outputpath = project.outputpath;
if iscell(filename)
    frames  = length(filename);
else
    frames = 1;
    filename = cellstr(filename);
end
set(gcf,'Pointer','watch');
for frame = 1 : frames
    sourcefile = fullfile(outputpath, [char(filename{frame}),'.mat']);
    destfile = fullfile(pathname,'snapshot',[char(filename{frame}),'.mat']);
    copyfile(sourcefile,destfile);
end
cd(currentpath);
CTsaveProject(handles,pathname,snapfilename);
set(gcf,'Pointer','arrow');
