function varargout = celldetectmodgui(varargin)
% CELLDETECTMODGUI M-file for celldetectmodgui.fig
%      CELLDETECTMODGUI, by itself, creates a new CELLDETECTMODGUI or raises the existing
%      singleton*.
%
%      H = CELLDETECTMODGUI returns the handle to a new CELLDETECTMODGUI or the handle to
%      the existing singleton*.
%
%      CELLDETECTMODGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLDETECTMODGUI.M with the given input arguments.
%
%      CELLDETECTMODGUI('Property','Value',...) creates a new CELLDETECTMODGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before celldetectgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to celldetectmodgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help celldetectmodgui

% Last Modified by GUIDE v2.5 22-Jan-2009 17:03:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @celldetectmodgui_OpeningFcn, ...
                   'gui_OutputFcn',  @celldetectmodgui_OutputFcn, ...
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


% --- Executes just before celldetectmodgui is made visible.
function celldetectmodgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to celldetectmodgui (see VARARGIN)
handles.Img_in=varargin{1}';
handles.historycount=0;
handles.histval=handles.historycount;
handles.corr_threshold=0.7;
set(handles.corrthresh,'String',num2str(handles.corr_threshold));
handles.reg_maxdist=0.00001;
set(handles.regmaxdist,'String',num2str(handles.reg_maxdist));

handles.strelval=1;
set(handles.strel,'String',num2str(handles.strelval));
handles.xdstval=1;
set(handles.xdst,'String',num2str(handles.xdstval));
handles.ydstval=1;
set(handles.ydst,'String',num2str(handles.ydstval));
handles.filtradval=2;
set(handles.filterradius,'String',num2str(handles.filtradval));
handles=estimate(handles);  
try
    handles.L=varargin{2};
catch
end
handles.L_prev=handles.L;
handles=scale(handles);
handles=refresh(handles);

% Choose default command line output for celldetectmodgui
handles.output = handles.L_out;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes celldetectmodgui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = celldetectmodgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);figure;

% --- Executes on slider movement.
function ar_lowth_Callback(hObject, eventdata, handles)
% hObject    handle to ar_lowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.ar_lowthtext,'String',num2str(round(get(handles.ar_lowth,'Value'))));
handles=refresh(handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ar_lowth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ar_lowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ecc_Callback(hObject, eventdata, handles)
% hObject    handle to ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.ecctext,'String',num2str(get(handles.ecc,'Value')));
handles=refresh(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ecc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ecc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ext_Callback(hObject, eventdata, handles)
% hObject    handle to ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.exttext,'String',num2str(get(handles.ext,'Value')));
handles=refresh(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ar_highth_Callback(hObject, eventdata, handles)
% hObject    handle to ar_highth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.ar_highthtext,'String',num2str(round(get(handles.ar_highth,'Value'))));
handles=refresh(handles);
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function ar_highth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ar_highth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in DONE.
function DONE_Callback(hObject, eventdata, handles)
% hObject    handle to DONE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = handles.L_out;
guidata(hObject, handles);
uiresume(handles.figure1)

function ecctext_Callback(hObject, eventdata, handles)
% hObject    handle to ecctext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ecctext as text
%        str2double(get(hObject,'String')) returns contents of ecctext as a double
val = str2double(get(hObject,'String'));
% Determine whether val is a number between 0 and 1.
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.ecc,'Min') && ...
   val <= get(handles.ecc,'Max')
   set(handles.ecc,'Value',val);
end
handles=refresh(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ecctext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ecctext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function exttext_Callback(hObject, eventdata, handles)
% hObject    handle to exttext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exttext as text
%        str2double(get(hObject,'String')) returns contents of exttext as a double
val = str2double(get(hObject,'String'));
% Determine whether val is a number between 0 and 1.
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.ext,'Min') && ...
   val <= get(handles.ext,'Max')
   set(handles.ext,'Value',val);
end
handles=refresh(handles);
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function exttext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exttext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ar_lowthtext_Callback(hObject, eventdata, handles)
% hObject    handle to ar_lowthtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ar_lowthtext as text
%        str2double(get(hObject,'String')) returns contents of ar_lowthtext as a double
val = str2double(get(hObject,'String'));
% Determine whether val is a number between 0 and 1.
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.ar_lowth,'Min') && ...
   val <= get(handles.ar_lowth,'Max')
   set(handles.ar_lowth,'Value',val);
end
handles=refresh(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ar_lowthtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ar_lowthtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ar_highthtext_Callback(hObject, eventdata, handles)
% hObject    handle to ar_highthtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ar_highthtext as text
%        str2double(get(hObject,'String')) returns contents of ar_highthtext as a double
val = str2double(get(hObject,'String'));
% Determine whether val is a number between 0 and 1.
if isnumeric(val) && length(val)==1 && ...
   val >= get(handles.ar_highth,'Min') && ...
   val <= get(handles.ar_highth,'Max')
   set(handles.ar_highth,'Value',val);
end
handles=refresh(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ar_highthtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ar_highthtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addroi.
function addroi_Callback(hObject, eventdata, handles)
% hObject    handle to addroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roifig=figure('Position',get(0,'ScreenSize'));
handles=refresh(handles);
while ishandle(roifig)
    try
        h=imellipse(gca);
        wait(h);
        roi=createMask(h);
        handles.L=handles.L|roi;
        handles=scale(handles);
        handles=refresh(handles);
        guidata(hObject, handles);
    catch
    end
end
handles=refresh(handles);

% --- Executes on button press in removeroi.
function removeroi_Callback(hObject, eventdata, handles)
% hObject    handle to removeroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roifig=figure('Position',get(0,'ScreenSize'));
handles=refresh(handles);
while ishandle(roifig)
    try
        h=imrect(gca);    
        roi=createMask(h);
        delete(h);
        handles.L=handles.L&~roi;
        handles=scale(handles);
        handles=refresh(handles);
        guidata(hObject, handles);
    catch
    end
end
handles=refresh(handles);


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=estimate(handles);
handles=scale(handles);
handles=refresh(handles);
guidata(hObject, handles);


% --- Executes on button press in dilate.
function dilate_Callback(hObject, eventdata, handles)
% hObject    handle to dilate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.L=imdilate(handles.L,strel('disk',handles.strelval));
handles=scale(handles);
handles=refresh(handles);
guidata(hObject, handles);


% --- Executes on button press in erode.
function erode_Callback(hObject, eventdata, handles)
% hObject    handle to erode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.L=imerode(handles.L,strel('disk',handles.strelval));
handles.cr=handles.cr-handles.strelval;
handles=scale(handles);
handles=refresh(handles);
guidata(hObject, handles);

% --- Executes on button press in rescale.
function rescale_Callback(hObject, eventdata, handles)
% hObject    handle to rescale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=scale(handles);
handles=refresh(handles);
guidata(hObject, handles);


% --- Executes on button press in addmagicroi.
function addmagicroi_Callback(hObject, eventdata, handles)
% hObject    handle to addmagicroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roifig=figure('Position',get(0,'ScreenSize'));
handles=refresh(handles);
while ishandle(roifig)
    
 try
p=impoint(gca);
pos=round(getPosition(p));
delete(p);
if logical(handles.Img_corr(pos(2),pos(1)))
J=regiongrowing(handles.bg_rec,pos(2),pos(1),handles.reg_maxdist);
handles.L=handles.L|logical(J);
handles=scale(handles);
handles=refresh(handles);
guidata(hObject, handles);

end
 catch
 end
end
handles=refresh(handles);

% --- Executes on button press in breakcells.
function breakcells_Callback(hObject, eventdata, handles)
% hObject    handle to breakcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roifig=figure('Position',get(0,'ScreenSize'));
handles=refresh(handles);
while ishandle(roifig)
    try
        p=imrect(gca);
        pos=round(getPosition(p));
        msk=createMask(p);
        L_out=breakcells(handles.L.*msk);
        handles.L(msk)=L_out(msk);
        delete(p);
        handles=scale(handles);
        handles=refresh(handles);
        guidata(hObject, handles);
    catch
    end
end
handles=refresh(handles);


function corrthresh_Callback(hObject, eventdata, handles)
% hObject    handle to corrthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corrthresh as text
%        str2double(get(hObject,'String')) returns contents of corrthresh as a double
val = str2double(get(hObject,'String'));
set(hObject,'String',num2str(val));
handles.corr_threshold=val;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function corrthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corrthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function regmaxdist_Callback(hObject, eventdata, handles)
% hObject    handle to regmaxdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of regmaxdist as text
%        str2double(get(hObject,'String')) returns contents of regmaxdist as a double
val = str2double(get(hObject,'String'))
set(hObject,'String',num2str(val));
handles.reg_maxdist=val;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function regmaxdist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regmaxdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function strel_Callback(hObject, eventdata, handles)
% hObject    handle to strel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strel as text
%        str2double(get(hObject,'String')) returns contents of strel as a double
val = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val));
handles.strelval=val;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function strel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles=refresh(handles)

handles.L_out=regionpropscheck(bwlabel(handles.L),get(handles.ar_lowth,'Value'),...
    get(handles.ar_highth,'Value'),get(handles.ecc,'Value'),...
    get(handles.ext,'Value'),[]);

imshow(handles.Img_ffcorr.*~bwperim(handles.L_out),...
    [min(min(handles.Img_ffcorr)) max(max(handles.Img_ffcorr))],...
    'InitialMagnification','fit');
handles.historycount=handles.historycount+1;
handles.histval=handles.historycount;
handles.Lhistory(:,:,handles.historycount)=handles.L_out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles=scale(handles)
props=regionprops(bwlabel(handles.L),'Area','Eccentricity','Extent');
for k=1:length(props)
    ar(k)=props(k).Area;
    ecc(k)=props(k).Eccentricity;
    ext(k)=props(k).Extent;
end
min_ar=min(ar);
max_ar=max(ar);
min_ecc=min(ecc);
max_ecc=max(ecc);
min_ext=min(ext);
max_ext=max(ext);
set(handles.ar_highth,'Max',max_ar,'Min',min_ar,'Value',max_ar);
set(handles.ar_highthtext,'String',num2str(max_ar));
set(handles.ar_lowth,'Max',max_ar,'Min',min_ar,'Value',min_ar);
set(handles.ar_lowthtext,'String',num2str(min_ar));
set(handles.ecc,'Max',1,'Min',0,'Value',1);
set(handles.ecctext,'String',num2str(1));
set(handles.ext,'Max',1,'Min',0,'Value',0);
set(handles.exttext,'String',num2str(0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles=estimate(handles)

handles.Img_ffcorr=ffcorr(handles.Img_in,50,0);% flat-field correction
[handles.Img_corr,handles.cr]=imgcorr(handles.Img_ffcorr);% image cross correlation
handles.Img_corr(handles.Img_corr<handles.corr_threshold)=0;
handles.Img_filt=ownmedfilt2tst(handles.Img_ffcorr,handles.filtradval);
[handles.bg_subtract, handles.bg_rec]=subtractBG(handles.Img_filt,handles.cr,0);
handles.seedstats=regionprops(bwlabel(imregionalmax(handles.bg_rec.*logical(handles.Img_corr))),'PixelList');
L=logical(zeros(size(handles.Img_in)));

for i=1:length(handles.seedstats)
   x=handles.seedstats(i).PixelList(1,2);
   y=handles.seedstats(i).PixelList(1,1);
   J=regiongrowing(handles.bg_rec,x,y,handles.reg_maxdist);
   L=L|logical(J);
end

handles.L=L;


% --- Executes on button press in excl_edges.
function excl_edges_Callback(hObject, eventdata, handles)
% hObject    handle to excl_edges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.L=(exclEdges(handles.L',handles.xdstval,handles.ydstval))';
handles=scale(handles);
        handles=refresh(handles);
        guidata(hObject, handles);

function xdst_Callback(hObject, eventdata, handles)
% hObject    handle to xdst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xdst as text
%        str2double(get(hObject,'String')) returns contents of xdst as a double
val = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val));
handles.xdstval=val;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xdst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xdst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ydst_Callback(hObject, eventdata, handles)
% hObject    handle to ydst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ydst as text
%        str2double(get(hObject,'String')) returns contents of ydst as a double
val = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val));
handles.ydstval=val;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ydst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ydst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function filterradius_Callback(hObject, eventdata, handles)
% hObject    handle to filterradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterradius as text
%        str2double(get(hObject,'String')) returns contents of filterradius as a double
val = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val));
handles.filtradval=val;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function filterradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in backbutton.
function backbutton_Callback(hObject, eventdata, handles)
% hObject    handle to backbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.histval>1
handles.histval=handles.histval-1;
handles.L_out=squeeze(handles.Lhistory(:,:,handles.histval));
imshow(handles.Img_ffcorr.*~bwperim(handles.L_out),...
    [min(min(handles.Img_ffcorr)) max(max(handles.Img_ffcorr))],...
    'InitialMagnification','fit');
guidata(hObject, handles);
end

% --- Executes on button press in forwardbutton.
function forwardbutton_Callback(hObject, eventdata, handles)
% hObject    handle to forwardbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.histval<size(handles.Lhistory,3)
handles.histval=handles.histval+1;
handles.L_out=squeeze(handles.Lhistory(:,:,handles.histval));
imshow(handles.Img_ffcorr.*~bwperim(handles.L_out),...
    [min(min(handles.Img_ffcorr)) max(max(handles.Img_ffcorr))],...
    'InitialMagnification','fit');
guidata(hObject, handles);
end

% --- Executes on button press in activatebutton.
function activatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to activatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.L=squeeze(handles.Lhistory(:,:,handles.histval));
handles=scale(handles);
        handles=refresh(handles);
        guidata(hObject, handles);