function [im mask] = CTreadImage(proj,frame,viewindex)
filename = proj.filename;
if iscell(filename)
    frames  = length(filename);
else
    frames = 1;
    filename = cellstr(filename);
end
outputpath = proj.outputpath;
currentpath = pwd;
cd(outputpath);
matfile = [char(filename{frame}),'.mat'];
load(matfile);
if viewindex <1
    viewindex = length(views);
end
data = views{viewindex};
im = data.im;
mask = data.mask;
cd(currentpath);            