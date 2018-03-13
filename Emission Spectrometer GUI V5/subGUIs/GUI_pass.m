function varargout = GUI_pass(varargin)
% GUI_PASS MATLAB code for GUI_pass.fig
%      GUI_PASS, by itself, creates a new GUI_PASS or raises the existing
%      singleton*.
%
%      H = GUI_PASS returns the handle to a new GUI_PASS or the handle to
%      the existing singleton*.
%
%      GUI_PASS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PASS.M with the given input arguments.
%
%      GUI_PASS('Property','Value',...) creates a new GUI_PASS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_pass_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_pass_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_pass

% Last Modified by GUIDE v2.5 09-Mar-2018 15:44:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_pass_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_pass_OutputFcn, ...
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


% --- Executes just before GUI_pass is made visible.
function GUI_pass_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_pass (see VARARGIN)

% Choose default command line output for GUI_pass
handles.output = 2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_pass wait for user response (see UIRESUME)
uiwait(handles.gui_pass_fig);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_pass_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.gui_pass_fig);


% --- Executes on button press in yes_btn.
function yes_btn_Callback(hObject, eventdata, handles)
% hObject    handle to yes_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 1;
uiresume(handles.gui_pass_fig);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in no_button.
function no_button_Callback(hObject, eventdata, handles)
% hObject    handle to no_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 2;
uiresume(handles.gui_pass_fig);
% Update handles structure
guidata(hObject, handles);




% --- Executes when user attempts to close gui_pass_fig.
function gui_pass_fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to gui_pass_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% The figure can be deleted now
delete(handles.gui_pass_fig);


% --- Executes on button press in no_close.
function no_close_Callback(hObject, eventdata, handles)
% hObject    handle to no_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = 3;
uiresume(handles.gui_pass_fig);
% Update handles structure
guidata(hObject, handles);
