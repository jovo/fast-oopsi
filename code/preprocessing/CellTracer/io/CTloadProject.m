function project = CTloadProject(h,pathname,filename)
handles = guidata(h);
currentpath = pwd;
cd(pathname);
load(filename,'-mat');
cd(currentpath);

project = proj;

% handles.viewtype = proj.viewtype;
% handles.viewindex = proj.viewindex;
% handles.project = proj.project;
% handles.io = proj.io;
% handles.link = proj.link;
% handles.undo = proj.undo; 
% 
% [im mask] = CTreadImage(proj.project,1,1);
% ax1 = subplot('position',[0.005, 0.20, 0.375 0.9]); 
% ax2 = subplot('position',[0.385, 0.20, 0.375 0.9]);
% ax1 = getAxes(gcf,1);subplot(ax1);
% imshow(im,[]); 
% ax2 = getAxes(gcf, 2);subplot(ax2);
% imshow(uint8(mask) * 255,[]);
% linkaxes([ax1 ax2]);
