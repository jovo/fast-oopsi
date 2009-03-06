function  setCurrentImage(h,im,mask,frame)
handles = guidata(h);
proj = handles.project;
if isempty(frame) || frame < 1
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
end
viewtype = handles.viewtype; %get which type of data to show 
if strcmp(viewtype,'input') %input
    CTsetImage(proj,im,mask,frame,handles.viewindex,1);
elseif strcmp(viewtype,'segmentation') %segmentation
    CTsetImage(proj,im,mask,frame,handles.viewindex,2);
elseif strcmp(viewtype,'tracking')
    CTsetImage(proj,[],im,frame,handles.viewindex,3);
    if frame < handles.io.frames
        CTsetImage(proj,[],mask,frame+1,handles.viewindex,3);
    end
end
guidata(h,handles);

function CTsetImage(proj,im,mask,frame,viewindex,viewtype)
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
if viewindex <1
    viewindex = length(views);
end
data = views{viewindex};
if viewtype ==1
    data.im = im;
    if ~isempty(mask)
        data.mask = mask>0;
    else
        data.mask = [];
    end
elseif viewtype ==2
    if ~isempty(mask)
        temp = data.mask;
        temp(:,:,3) = mask>0;
        data.mask = temp;
    end
elseif viewtype == 3
    if ~isempty(im)
        data.im = im;
    end
    if ~isempty(mask)
        if length(size(mask)) == 3
            data.mask = mask;
        else
            temp = data.mask;
            temp(:,:,3) = mask>0;
            data.mask = temp;
        end
    end
end
views{viewindex} = data;
save(matfile, 'views');
cd(currentpath);          