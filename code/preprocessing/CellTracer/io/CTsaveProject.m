function CTsaveProject(handles,pathname,filename)
proj.viewtype = handles.viewtype;
proj.viewindex = handles.viewindex;
proj.project = handles.project;
proj.io = handles.io;
proj.link = handles.link;
proj.undo = handles.undo;

currentpath = pwd;
cd(pathname);
save(filename, 'proj');
cd(currentpath);