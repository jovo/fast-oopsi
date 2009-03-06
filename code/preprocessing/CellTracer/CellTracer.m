function varargout = CellTracer(varargin)
% CELLTRACER M-file for CellTracer.fig
%      CELLTRACER, by itself, creates a new CELLTRACER or raises the existing
%      singleton*.
%
%      H = CELLTRACER returns the handle to a new CELLTRACER or the handle to
%      the existing singleton*.
%
%      CELLTRACER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLTRACER.M with the given input arguments.
%
%      CELLTRACER('Property','Value',...) creates a new CELLTRACER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellTracer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellTracer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Edit menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellTracer

% Last Modified by GUIDE v2.5 10-Jun-2008 09:44:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellTracer_OpeningFcn, ...
                   'gui_OutputFcn',  @CellTracer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before CellTracer is made visible.
function CellTracer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellTracer (see VARARGIN)

% Choose default command line output for CellTracer
handles.output = hObject;
addpath('gui','io','resource','seg','tracking','lines');
handles.viewtype = 'input';   %1 for input, 2 for segmentaion, 3 for tracking
handles.viewindex = 1;
handles.undo.flag = 0;   % no undo at the begining
handles.io.frames = 0;
handles.project = [];
guidata(hObject, handles);
addToolBar();
updateMenuStatus();
% Update handles structure

% UIWAIT makes CellTracer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CellTracer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function NewMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to NewMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.jpg;*.tif;*.tiff;*.png;*.bmp','Image Files';
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file',...
   'MultiSelect', 'on');
% Update handles structure
if isequal(filename,0)
   set(handles.StatusBar,'String', 'Reading input image files....Canceled.');
else
    handles.project.tag = 'project1';
    handles.project.filename = filename;
    handles.project.pathname = pathname;
    handles.project.outputpath = [pathname,'output'];
    set(handles.StatusBar,'String', 'Reading input image files....');
    
    [handles.io.frames statustext]= getImages(filename,pathname,handles.project.outputpath,hObject);
    handles.link = {};
    set(handles.StatusBar,'String', statustext);
    if handles.io.frames>0
        set(handles.FrameIndex, 'String',num2str(handles.io.frames));
        updateButtons(handles.io.frames,handles.io.frames,handles);
    else
        set(handles.FrameIndex, 'String','');
    end
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.mat','Project File'}, ...
   'Load an existing project',...
   'MultiSelect', 'off');
if isequal(filename,0)
   set(handles.StatusBar,'String', 'Reading input project file....Canceled.');
else
    proj = CTloadProject(hObject,pathname,filename);
    %give user a chance to correct file path problem
    button = questdlg(['Is image file path ',proj.project.pathname, ' correct?'],'Image File Path Confirmation');
    if ~strcmp(button,'Yes')
        proj.project.pathname = [uigetdir,filesep];
        proj.project.outputpath = [proj.project.pathname, 'output'];
    end
    
    %update project infomation
    handles.viewtype = 'input';
    handles.viewindex = proj.viewindex;
    handles.project = proj.project;
    handles.io = proj.io;
    handles.link = proj.link;
    handles.undo = proj.undo; 
    CTreloadImage(handles.project,hObject);
    if handles.io.frames>0
        set(handles.FrameIndex, 'String',num2str(1));
        updateButtons(1,handles.io.frames,handles);
        %showFrame(hObject,1);
        set(handles.StatusBar,'String', ['Loading project done. Number of frames:' num2str( handles.io.frames)]);
    else
        set(handles.FrameIndex, 'String','');
        set(handles.StatusBar,'String', 'Loading project failed.');
    end
    guidata(hObject,handles);
    setView(hObject,'input',1);
end

% --------------------------------------------------------------------
function menuExportMask_Callback(hObject, eventdata, handles)
% hObject    handle to menuExportMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for frame = 1:handles.io.frames
    CTexportImage(handles.project,frame,2);
end
set(handles.StatusBar,'String', 'Exporting segmentation results done.');

% --------------------------------------------------------------------
function menuExportImage_Callback(hObject, eventdata, handles)
% hObject    handle to menuExportImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('This operation will replace input images with the modified images. Do you want to continue?','Confirmation');
if strcmp(button,'Yes')
    for frame = 1:handles.io.frames
    CTexportImage(handles.project,frame,1);
end
set(handles.StatusBar,'String', 'Exporting images done.');
end


% --------------------------------------------------------------------
function ExitMenuItem_Callback(hObject, eventdata, handles)
CTsaveProject(handles,'.','last.mat');
close
 
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.mat';},'Save Project as');
CTsaveProject(handles,pathname,filename);

% --------------------------------------------------------------------
function menuRestoreSnapshot_Callback(hObject, eventdata, handles)
% hObject    handle to menuRestoreSnapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.sns','Project Snap File'}, ...
   'Load an existing project snapshot',...
   'MultiSelect', 'off');
if isequal(filename,0)
   set(handles.StatusBar,'String', 'Reading input project file....Canceled.');
else
    CTrestoreFromSnapshot(hObject,pathname,filename);
    handles = guidata(hObject);
    handles.viewtype = 'input';
    guidata(hObject,handles);
    if handles.io.frames>0
        set(handles.FrameIndex, 'String',num2str(1));
        updateButtons(1,handles.io.frames,handles);
        %showFrame(hObject,-1);
        set(handles.StatusBar,'String', ['Restore project snapshot done. Number of frames:' num2str( handles.io.frames)]);
    else
        set(handles.FrameIndex, 'String','');
        set(handles.StatusBar,'String', 'Loading project snapshot failed.');
    end
    setView(hObject,'input',1);
end
% --------------------------------------------------------------------
function menuSnapshot_Callback(hObject, eventdata, handles)
% hObject    handle to menuSnapshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.sns';},'Save Project Snapshot as');
CTSnapshot(hObject,pathname,filename);

% --------------------------------------------------------------------
function Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in FirstButton.
function FirstButton_Callback(hObject, eventdata, handles)
% hObject    handle to FirstButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame= 1;
updateButtons(frame, handles.io.frames,handles);
set(handles.FrameIndex,'String',num2str(frame));
showFrame(hObject,frame);

% --- Executes on button press in PreviousButton.
function PreviousButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
frame = frame - 1;
updateButtons(frame, handles.io.frames,handles);
set(handles.FrameIndex,'String',num2str(frame));
showFrame(hObject,frame);

% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
frame = frame + 1;
updateButtons(frame, handles.io.frames,handles);
set(handles.FrameIndex,'String',num2str(frame));
showFrame(hObject,frame);

% --- Executes on button press in LastButton.
function LastButton_Callback(hObject, eventdata, handles)
% hObject    handle to LastButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame= handles.io.frames;
updateButtons(frame, handles.io.frames,handles);
set(handles.FrameIndex,'String',num2str(frame));
showFrame(hObject,frame);

% --- Executes on button press in RedrewButton.
function RedrewButton_Callback(hObject, eventdata, handles)
% hObject    handle to RedrewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
showFrame(hObject,frame);


% --- Executes on button press in PlayButton.
function PlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frames = handles.io.frames;
for frame = 1: frames
    updateButtons(frame, frames,handles);
    set(handles.FrameIndex,'String',num2str(frame));
    showFrame(hObject,frame);
    pause(0.25);
end


% --- Executes on button press in FirstCell.
function FirstCell_Callback(hObject, eventdata, handles)
% hObject    handle to FirstCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cell= 1; maxcell = getMaxCellIndex(hObject,[]);
updateCellButtons(cell, maxcell,handles);
set(handles.CellIndex,'String',num2str(cell));
showCell(hObject,cell);


% --- Executes on button press in PreviousCell.
function PreviousCell_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cell= get(handles.CellIndex,'String'); cell = str2num(cell);
cell = cell - 1;maxcell = getMaxCellIndex(hObject,[]);
updateCellButtons(cell, maxcell,handles);
set(handles.CellIndex,'String',num2str(cell));
showCell(hObject,cell);
guidata(hObject, handles);


% --- Executes on button press in NextCell.
function NextCell_Callback(hObject, eventdata, handles)
% hObject    handle to NextCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cell= get(handles.CellIndex,'String'); cell = str2num(cell);
cell = cell + 1;maxcell = getMaxCellIndex(hObject,[]);
updateCellButtons(cell, maxcell,handles);
set(handles.CellIndex,'String',num2str(cell));
showCell(hObject,cell);
guidata(hObject, handles);

% --- Executes on button press in LastCell.
function LastCell_Callback(hObject, eventdata, handles)
% hObject    handle to LastCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxcell = getMaxCellIndex(hObject,[]); cell = maxcell;
updateCellButtons(cell, maxcell,handles);
set(handles.CellIndex,'String',num2str(cell));
showCell(hObject,cell);
guidata(hObject, handles);

% --- Executes on button press in PlayCell.
function PlayCell_Callback(hObject, eventdata, handles)
% hObject    handle to PlayCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxcell = getMaxCellIndex(hObject,[]);
for cell = 1: maxcell
    updateCellButtons(cell, maxcell,handles);
    set(handles.CellIndex,'String',num2str(cell));
    showCell(hObject,cell);
    pause(0.25);
end
guidata(hObject, handles);

% --- Executes on button press in RedrawCell.
function RedrawCell_Callback(hObject, eventdata, handles)
% hObject    handle to RedrawCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cell= get(handles.CellIndex,'String'); cell = str2num(cell);
maxcell = getMaxCellIndex(hObject,[]);
if cell > maxcell
    cell = maxcell;
    updateCellButtons(cell, maxcell,handles);
    set(handles.CellIndex,'String',num2str(cell));
end
showCell(hObject,cell);
guidata(hObject, handles);

function CellIndex_Callback(hObject, eventdata, handles)
% hObject    handle to CellIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellIndex as text
%        str2double(get(hObject,'String')) returns contents of CellIndex as a double
cell= get(handles.CellIndex,'String');
maxcell = getMaxCellIndex(hObject,[]);
%convert from string to number if possible, otherwise returns empty
try
    cell = str2num(cell);
    if cell > maxcell || cell < 1
        cell = maxcell;
        set(handles.CellIndex,'String',num2str(cell));
        set(handles.StatusBar,'String', ['Invalid input: should be a number between 1 and ' num2str(cell)]);
    end
    showCell(hObject,cell);
    updateCellButtons(cell,maxcell,handles);
catch
    cell = maxcell;
    set(handles.CellIndex,'String',num2str(cell));
    showCell(hObject,cell);
    updateButtons(cell,maxcell,handles);
    set(handles.StatusBar,'String', 'Input invalid');
end
guidata(hObject,handles);

function FrameIndex_Callback(hObject, eventdata, handles)
% hObject    handle to FrameIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameIndex as text
%        str2double(get(hObject,'String')) returns contents of FrameIndex as a double
%get the string for the editText component
frame= get(handles.FrameIndex,'String');
%convert from string to number if possible, otherwise returns empty
try
    frame = str2num(frame);
    if frame > handles.io.frames || frame < 1
        frame = handles.io.frames;
        set(handles.FrameIndex,'String',num2str(frame));
        set(handles.StatusBar,'String', ['Invalid input: should be a number between 1 and ' num2str(frame)]);
    end
    showFrame(hObject,frame);
    updateButtons(frame,handles.io.frames,handles);
catch
    frame = handles.io.frames;
    set(handles.FrameIndex,'String',num2str(frame));
    showFrame(hObject,frame);
    updateButtons(frame,frame,handles);
    set(handles.StatusBar,'String', 'Input invalid');
end

% --------------------------------------------------------------------
function menuDoubleResolution_Callback(hObject, eventdata, handles)
% hObject    handle to menuDoubleResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Double Resolution?','Confirmation');
%if handles.seg.lastversion == 0
    if strcmp(button,'Yes')
        frames = handles.io.frames;
        ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
        handles.undo.flag = 1;
        %handles.undo.data = handles.im{handles.io.index};
        handles.undo.type = 'input';
        for frame =1 : frames
            [temp tempmask]=getCurrentImage(hObject,frame);
            im = ExpandFrame(temp);
            if ~isempty(tempmask)
                mask = [];
                mask(:,:,1) = im<0;
                mask(:,:,2) = im<0;
                mask(:,:,3) = im<0;
            end
            setCurrentImage(hObject,im,mask,frame);
            purgeCurrentImage(hObject,frame,1);%%% to set to the real viewindex later
        end
        guidata(hObject,handles);
        showFrame(hObject,-1); %without rescale
        updateMenuStatus();
    end
%end

% --------------------------------------------------------------------
function menutrackingresetall_Callback(hObject, eventdata, handles)
% hObject    handle to menutrackingresetall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Reset all tracking results?','Confirmation');
if strcmp(button,'Yes')
    frames = handles.io.frames;
    for frame =1 : frames
        %[clink ccolor oldalignment] = getCurrentLink(hObject,[]);
        setCurrentLink(hObject,[],[],[],frame);
    end
    showFrame(hObject,[]); %without rescale
    updateMenuStatus();
end


% --------------------------------------------------------------------
function menutrackingresetlineage_Callback(hObject, eventdata, handles)
% hObject    handle to menutrackingresetlineage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Reset all tracking results?','Confirmation');

if strcmp(button,'Yes')
    frames = handles.io.frames;
    for frame =1 : frames
        [clink ccolor oldalignment] = getCurrentLink(hObject,[]);
        setCurrentLink(hObject,[],[],oldalignment,frame);
    end
    showFrame(hObject,[]); 
    updateMenuStatus();
end


% --------------------------------------------------------------------
function menuCrop_Callback(hObject, eventdata, handles)
% hObject    handle to menuCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Crop images to current view?','Crop confirmation');
%if handles.seg.lastversion == 0
    if strcmp(button,'Yes')
        frames = handles.io.frames;
        [temp tempmask]=getCurrentImage(hObject,1);
        ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
        [m n] = size(temp);
        v = uint16(v);v(v<1) = 1; 
        if v(2) > n v(2) = n; end
        if v(4) > m v(4) = m; end
        handles.undo.flag = 1;
        %handles.undo.data = handles.im{handles.io.index};
        handles.undo.type = 'input';
        
        for frame =1 : frames
            [temp tempmask]=getCurrentImage(hObject,frame);
            im = temp(v(3):v(4),v(1):v(2));
            mask = tempmask(v(3):v(4),v(1):v(2),:);
            setCurrentImage(hObject,im,mask,frame);
            purgeCurrentImage(hObject,frame,1);%%% to set to the real viewindex later
        end
        guidata(hObject,handles);
        imshow(temp);
        showFrame(hObject,-1);
        updateMenuStatus();
    end
%end

% --------------------------------------------------------------------
function menuInvert_Callback(hObject, eventdata, handles)
% hObject    handle to menuInvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Invert all images?','Invert confirmation');
%if handles.seg.lastversion == 0
    if strcmp(button,'Yes')
        frames = handles.io.frames;
        ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
        handles.undo.flag = 1;
        %handles.undo.data = handles.im{handles.io.index};
        handles.undo.type = 'input';
        
        for frame =1 : frames
            [temp tempmask]=getCurrentImage(hObject,frame);
            im = 255 - temp; im(im==0) = 1;
            mask = [];
            setCurrentImage(hObject,im,mask,frame);
            purgeCurrentImage(hObject,frame,1);%%% to set to the real viewindex later
        end
        guidata(hObject,handles);
        imshow(temp);
        showFrame(hObject,[]);
        updateMenuStatus();
    end
%end


function cmenu_removemaskedobject(hObject, eventdata, handles)
userdata = get(hObject,'UserData');
pts = userdata{1};
currentview = userdata{2};
[im mask] = getCurrentImage(hObject,[]);
if currentview == 1 % can only deal with right panel
    if length(size(mask)) == 3
        [bwmask1 nummask1] = bwlabel(mask(:,:,1),4); 
        [bwmask2 nummask2] = bwlabel(mask(:,:,2),4);
        [bwmask3 nummask3] = bwlabel(mask(:,:,3),4);
        nummask = nummask1 + nummask2 + nummask3;
        bwmask(:,:,1) = bwmask1;
        bwmask(:,:,2) = bwmask2;
        bwmask(:,:,3) = bwmask3;
    else
        bwmask = mask;
        [bwmask nummask] = bwlabel(bwmask,4); 
    end
    
    pts =pts{2}; [pm pn] = size(pts);
    for j = 1:pm
        if pts(j,1) > 0

            if length(size(bwmask)) == 3
                cellindex = bwmask(pts(j,3),pts(j,2),:);
            else
                cellindex = bwmask(pts(j,3),pts(j,2));
            end

            if length(size(mask)) == 3
                for k = 1:3
                    if cellindex(k) > 0
                        t = mask(:,:,k); t(bwmask(:,:,k) == cellindex(k)) = 0; mask(:,:,k) = t;
                    end

                end
            else
                mask(bwmask == cellindex) = 0;
            end
        end
    end
    setCurrentImage(hObject,im,mask,[]);
    v = axis; imshow(uint8(mask) * 255); axis(v);
elseif currentview == 2
    
    bwmask = mask(:,:,1) == 0 & mask(:,:,2) == 255 & mask(:,:,3) == 255; 
    [bwmask nummask] = bwlabel(~bwmask,4); 

    pts =pts{2}; [pm pn] = size(pts);
    for j = 1:pm
        if pts(j,1) > 0
            cellindex = bwmask(pts(j,3),pts(j,2));
            bwmask(bwmask==cellindex) = 0;
        end
    end
    setCurrentImage(hObject,im,bwmask,[]);
    showFrame(hObject,[]);
elseif currentview == 3
    handles = guidata(hObject);
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
    
    %deal with the left panel first
    bwmask = getCurrentLabels(hObject,frame);
    pts1 =pts{1}; [pm pn] = size(pts1);
    for j = 1:pm
        if pts1(j,1) > 0
            cellindex = bwmask(pts1(j,3),pts1(j,2));
            im(bwmask == cellindex) = 0;
            bwmask(bwmask==cellindex) = 0;
        end
    end
    setCurrentLabels(hObject,bwmask,frame);
    
    bwmask = getCurrentLabels(hObject,frame+1);
    if ~isempty(bwmask)
        pts2 =pts{2}; [pm pn] = size(pts2);
        for j = 1:pm
            if pts2(j,1) > 0
                cellindex = bwmask(pts2(j,3),pts2(j,2));
                mask(bwmask == cellindex) = 0;
                bwmask(bwmask==cellindex) = 0;
            end
        end
        setCurrentLabels(hObject,bwmask,frame+1);
    end
    setCurrentImage(hObject,im,mask,[]);
    showFrame(hObject,[]);
end

function cmenu_linkmaskedobject(hObject, eventdata, handles)
pts =get(hObject,'UserData');

[im mask] = getCurrentImage(hObject,[]);  %get curent images based on view
[clink ccolor alignment] = getCurrentLink(hObject,[]);  %get all verified links
[pm1 pn1] = size(pts{1});
[pm2 pn2] = size(pts{2});
if isempty(clink) %no links defined before, so create a link
    clink = zeros(pm1,pm2);
else 
    %unlink previous links 
    for i = 1:pm1
        if pts{1}(i,1) > 0
            clink(i,:) = 0;
        end
    end
    for i = 1:pm2
        if pts{2}(i,1) > 0
            clink(:,i) = 0;
        end
    end
end

%now add the new links
clink(pts{1}(:,1) > 0, pts{2}(:,1) > 0) = 1;

%now update the colors
ccolor = CTassignLinkColor(clink, ccolor);
setCurrentLink(hObject,clink,ccolor,alignment,[]);
redrawImages(im,mask,ccolor);

% --------------------------------------------------------------------
function View_Callback(hObject, eventdata, handles)
% hObject    handle to View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ViewInput_Callback(hObject, eventdata, handles)
% hObject    handle to ViewInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setView(hObject,'input',1);

% --------------------------------------------------------------------
function ViewSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to ViewSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setView(hObject,'segmentation',1);

% --------------------------------------------------------------------
function ViewTracking_Callback(hObject, eventdata, handles)
% hObject    handle to ViewTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setView(hObject,'tracking',1);
%showCell(hObject,1);


function setView(hObject,viewtype,viewindex)
% handles = guidata(hObject);
% handles.viewtype = viewtype;
% guidata(hObject,handles);
CTsetView(viewtype,viewindex);
updateMenuStatus();
showFrame(hObject,[]);



% --- Executes during object creation, after setting all properties.
function CellIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white menusegmentation on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menuUndo_Callback(hObject, eventdata, handles)
% hObject    handle to menuUndo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTUndo;

% --------------------------------------------------------------------
function menuSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to menuSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menubackground1_Callback(hObject, eventdata, handles)
% hObject    handle to menubackground1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('background1');

% --- Executes on button press in ButtonRun.
function ButtonRun_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
para = CTgetParameters;
method = get(hObject,'TooltipString');
switch method
    case {'background1'}
        CTrunInitialPartition(1,2);
    case {'background2'}
       CTrunInitialPartition(2,2);
    case {'background3'}
       CTrunInitialPartition(3,2);
    case {'border1'}
       CTrunInitialPartition(4,1);
    case {'border2'}
       CTrunInitialPartition(5,1);
    case {'border3'}
       CTrunInitialPartition(6,1);
    case {'border4'}
       CTrunInitialPartition(7,1);
    case {'cell1'}
       CTrunInitialPartition(10,3);
    case {'cell2'}
       CTrunInitialPartition(11,3);
    case {'globalalignment'}
        CTrunTracking(1);
    case {'celltracking'}
        CTrunTracking(2);
    otherwise
        msgbox('Unknown method.')
    
end


% --- Executes on button press in checkboxAllFrames.
function checkboxAllFrames_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAllFrames


% --- Executes on button press in checkboxDebug.
function checkboxDebug_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDebug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDebug




% % --------------------------------------------------------------------
% function MenuTracking_Callback(hObject, eventdata, handles)
% % hObject    handle to MenuTracking (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% %CTgetParameters;
% CTsetParameters('background1');


% --------------------------------------------------------------------
function menubackground2_Callback(hObject, eventdata, handles)
% hObject    handle to menubackground2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('background2');


% --------------------------------------------------------------------
function menubackground3_Callback(hObject, eventdata, handles)
% hObject    handle to menubackground3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('background3');

% --------------------------------------------------------------------
function menuborder1_Callback(hObject, eventdata, handles)
% hObject    handle to menuborder1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('border1');

% --------------------------------------------------------------------
function menuborder2_Callback(hObject, eventdata, handles)
% hObject    handle to menuborder2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('border2');

% --------------------------------------------------------------------
function menuborder3_Callback(hObject, eventdata, handles)
% hObject    handle to menuborder3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('border3');

% --------------------------------------------------------------------
function menuborder4_Callback(hObject, eventdata, handles)
% hObject    handle to menuborder4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('border4');
% --------------------------------------------------------------------
function menuCell1_Callback(hObject, eventdata, handles)
% hObject    handle to menuCell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('cell1');

% --------------------------------------------------------------------
function menucell2_Callback(hObject, eventdata, handles)
% hObject    handle to menucell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('cell2');

% --------------------------------------------------------------------
function menutracking1_Callback(hObject, eventdata, handles)
% hObject    handle to menutracking1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('globalalignment');

% --------------------------------------------------------------------
function menucelltracking_Callback(hObject, eventdata, handles)
% hObject    handle to menucelltracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTsetParameters('celltracking');

function editPrimary1_Callback(hObject, eventdata, handles)
% hObject    handle to editPrimary1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPrimary1 as text
%        str2double(get(hObject,'String')) returns contents of editPrimary1 as a double


% --- Executes during object creation, after setting all properties.
function editPrimary1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPrimary1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPrimary2_Callback(hObject, eventdata, handles)
% hObject    handle to editPrimary2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPrimary2 as text
%        str2double(get(hObject,'String')) returns contents of editPrimary2 as a double


% --- Executes during object creation, after setting all properties.
function editPrimary2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPrimary2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSecondary3_Callback(hObject, eventdata, handles)
% hObject    handle to editSecondary3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSecondary3 as text
%        str2double(get(hObject,'String')) returns contents of editSecondary3 as a double


% --- Executes during object creation, after setting all properties.
function editSecondary3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSecondary3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSecondary1_Callback(hObject, eventdata, handles)
% hObject    handle to editSecondary1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSecondary1 as text
%        str2double(get(hObject,'String')) returns contents of editSecondary1 as a double


% --- Executes during object creation, after setting all properties.
function editSecondary1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSecondary1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editSecondary2_Callback(hObject, eventdata, handles)
% hObject    handle to editSecondary2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSecondary2 as text
%        str2double(get(hObject,'String')) returns contents of editSecondary2 as a double


% --- Executes during object creation, after setting all properties.
function editSecondary2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSecondary2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSecondary4_Callback(hObject, eventdata, handles)
% hObject    handle to editSecondary4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSecondary4 as text
%        str2double(get(hObject,'String')) returns contents of editSecondary4 as a double


% --- Executes during object creation, after setting all properties.
function editSecondary4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSecondary4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSecondary5_Callback(hObject, eventdata, handles)
% hObject    handle to editSecondary5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSecondary5 as text
%        str2double(get(hObject,'String')) returns contents of editSecondary5 as a double


% --- Executes during object creation, after setting all properties.
function editSecondary5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSecondary5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxSequential.
function checkboxSequential_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSequential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSequential



function editPrimary3_Callback(hObject, eventdata, handles)
% hObject    handle to editPrimary3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPrimary3 as text
%        str2double(get(hObject,'String')) returns contents of editPrimary3 as a double

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function menuhelp_Callback(hObject, eventdata, handles)
% hObject    handle to menuhelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function MenuTracking_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function menuHelpAbout_Callback(hObject, eventdata, handles)
% hObject    handle to menuHelpAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rect = get(gcf,'Position');
cx = rect(1) + rect(3) / 2;
cy = rect(2) + rect(4) / 2;
h = CTabout;
rect =get(h,'Position');
rect(1) = cx-rect(3); rect(2) = cy + rect(4);
set(h,'Position',rect);

% --------------------------------------------------------------------
function menuHelpWeb_Callback(hObject, eventdata, handles)
% hObject    handle to menuHelpWeb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web 'www.stat.duke.edu/research/software/west/celltracer' -browser



% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --------------------------------------------------------------------
function munuPermentLabels_Callback(hObject, eventdata, handles)
% hObject    handle to munuPermentLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CTcacheLineageImages(hObject,[]);
