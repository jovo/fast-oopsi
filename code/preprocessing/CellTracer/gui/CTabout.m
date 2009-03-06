function varargout = CTabout(varargin)
% CTABOUT M-file for CTabout.fig
%      CTABOUT, by itself, creates a new CTABOUT or raises the existing
%      singleton*.
%
%      H = CTABOUT returns the handle to a new CTABOUT or the handle to
%      the existing singleton*.
%
%      CTABOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CTABOUT.M with the given input arguments.
%
%      CTABOUT('Property','Value',...) creates a new CTABOUT or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CTabout_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CTabout_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CTabout

% Last Modified by GUIDE v2.5 06-Jun-2008 08:00:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CTabout_OpeningFcn, ...
                   'gui_OutputFcn',  @CTabout_OutputFcn, ...
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

% --- Executes just before CTabout is made visible.
function CTabout_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CTabout (see VARARGIN)

% Choose default command line output for CTabout
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CTabout wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CTabout_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




