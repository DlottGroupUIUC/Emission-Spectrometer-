function varargout = binning_GUI(varargin)
% BINNING_GUI MATLAB code for binning_GUI.fig
%      BINNING_GUI, by itself, creates a new BINNING_GUI or raises the existing
%      singleton*.
%
%      H = BINNING_GUI returns the handle to a new BINNING_GUI or the handle to
%      the existing singleton*.
%
%      BINNING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BINNING_GUI.M with the given input arguments.
%
%      BINNING_GUI('Property','Value',...) creates a new BINNING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before binning_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to binning_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help binning_GUI

% Last Modified by GUIDE v2.5 12-Mar-2018 10:12:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @binning_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @binning_GUI_OutputFcn, ...
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


% --- Executes just before binning_GUI is made visible.
function binning_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to binning_GUI (see VARARGIN)

%load data from main GUI
h = findobj('Tag','Emission_spectrometer_GUI');
 
  % if exists (not empty)
 if ~isempty(h)
    % get handles and other user-defined data associated to Gui1
    main_GUI_data = guidata(h);

 end
 
%Set file index and initial file
file_idx=1;
file_list=main_GUI_data.file_list;
PMT_file=[file_list{file_idx},'.tdms'];

%Set table data
handles.binning_table.Data=main_GUI_data.central_table.Data(:,[1,4,5,6]);

%Set current plot
handles.selected_file_txt.String=main_GUI_data.file_list{file_idx};

%%%Initialize Plots%%%
%Generate data
time=logspace(-9,-5,1e3);
rad=zeros(1,1e3);
srad=zeros(32,1e3);

axes(handles.binning_axis);
handles.rad_sub=subplot('Position',[0.35 0.25 0.60 0.50]);handles.rad_plot=semilogx(time,rad,'b');
handles.srad_sub=subplot('Position',[0.35 0.75 0.60 0.20]);handles.srad_plot=semilogx(time,srad);

set(handles.srad_sub,'XTickLabel',[]);
ylabel(handles.rad_sub,'radiance (W/sr \cdot m^2)');
ylabel(handles.srad_sub,'spectral radiance (W/sr \cdot m^3)');
xlabel(handles.rad_sub,'time (s)');
xlim(handles.rad_sub,[1e-9 1e-5]);
linkaxes([handles.rad_sub handles.srad_sub],'x');


%Set filters to initialize
handles.emission_no_filter.Value=1;
handles.cal_no_filter.Value=1;

handles.emiss_filter_flag=0;
handles.cal_filter_flag=0;

%initialize supressed channels
handles.supress_channel.String='0';

%initialize target and threshold
handles.target_time.String='5';
handles.threshold_RMS.String='4';
handles.delay_flag=0;

%initialize close program decision
handles.decision = 3;

%Check for calibration
if size(main_GUI_data.main_data_table,2)==15 && isempty(main_GUI_data.main_data_table{1,15})~=1;
    handles.unit_conv=main_GUI_data.main_data_table{1,15};
end

% Update main GUI
handles.output = main_GUI_data;
handles.main_GUI_data=main_GUI_data;
handles.file_idx=file_idx;
handles.PMT_file=PMT_file;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes binning_GUI wait for user response (see UIRESUME)
% uiwait(handles.binning_GUI_bkg);


% --- Outputs from this function are returned to the command line.
function varargout = binning_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        varargout{1} = handles.output;

        

% --- Executes during object creation, after setting all properties.
function binning_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binning_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%initialize table
set(hObject,'Data',cell(1,4));
set(hObject, 'ColumnName', {'shot #','impact time (ns)','PMT','radiance' });
set(hObject, 'ColumnFormat',{'numeric','numeric','logical','logical',});


% --- Executes during object creation, after setting all properties.
function selected_file_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_file_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in p1_btn.
function p1_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p1_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;



handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+1;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+1;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+1;


table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);




% --- Executes on button press in p2_btn.
function p2_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p2_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+2;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+2;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+2;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);




% --- Executes on button press in p3_btn.
function p3_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p3_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+3;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+3;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+3;

table_single_Callback(hObject, eventdata, handles);
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);




% --- Executes on button press in p4_btn.
function p4_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p4_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+4;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+4;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+4;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);




% --- Executes on button press in p5_btn.
function p5_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p5_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+5;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+5;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+5;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in p10_btn.
function p10_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p10_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+10;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+10;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+10;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in p25_btn.
function p25_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p25_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+25;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+25;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+25;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in p50_btn.
function p50_btn_Callback(hObject, eventdata, handles)
% hObject    handle to p50_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}+50;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}+50;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}+50;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n50_btn.
function n50_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n50_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-50;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-50;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-50;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n25_btn.
function n25_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n25_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-25;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-25;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-25;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n10_btn.
function n10_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n10_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-10;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-10;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-10;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n5_btn.
function n5_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n5_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-5;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-5;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-5;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n4_btn.
function n4_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n4_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-4;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-4;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-4;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n3_btn.
function n3_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n3_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-3;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-3;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-3;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n2_btn.
function n2_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n2_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-2;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-2;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-2;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);


% --- Executes on button press in n1_btn.
function n1_btn_Callback(hObject, eventdata, handles)
% hObject    handle to n1_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;

handles.main_GUI_data.main_data_table{file_idx,4}=handles.main_GUI_data.main_data_table{file_idx,4}-1;
handles.binning_table.Data{file_idx,2}=handles.binning_table.Data{file_idx,2}-1;
handles.main_GUI_data.central_table.Data{file_idx,4}=handles.main_GUI_data.central_table.Data{file_idx,4}-1;

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);
guidata(hObject,handles);



function target_time_Callback(hObject, eventdata, handles)
% hObject    handle to target_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_time as text
%        str2double(get(hObject,'String')) returns contents of target_time as a double


% --- Executes during object creation, after setting all properties.
function target_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Table_header_Callback(hObject, eventdata, handles)
% hObject    handle to Table_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Peak_header_Callback(hObject, eventdata, handles)
% hObject    handle to Peak_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function rise_header_Callback(hObject, eventdata, handles)
% hObject    handle to rise_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function rise_single_Callback(hObject, eventdata, handles)
% hObject    handle to rise_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;PMT_p=handles.main_GUI_data.PMT_p;
file_list=handles.main_GUI_data.file_list;

%Generate PMT file name
PMT_file=[file_list{file_idx} '.tdms'];

%Grab target and threshold 
target= str2num(handles.target_time.String)*1e-9;
threshold= str2num(handles.threshold_RMS.String);

%Calculate flyer delay
rise_delay=peakdelay_rising_edge_v5(PMT_file,PMT_p,target,threshold);

%Update tables with delay
handles.main_GUI_data.main_data_table{file_idx,4}=round(rise_delay*1e9,0);
handles.binning_table.Data{file_idx,2}=round(rise_delay*1e9,0);
handles.main_GUI_data.central_table.Data{file_idx,4}=round(rise_delay*1e9,0);


handles.delay_flag=1;
handles.rise_delay=rise_delay;
handles.rise_idx=file_idx;
guidata(hObject,handles);

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);

handles.delay_flag=0;
guidata(hObject,handles);

% --------------------------------------------------------------------
function rise_all_Callback(hObject, eventdata, handles)
% hObject    handle to rise_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_p=handles.main_GUI_data.PMT_p;
file_list=handles.main_GUI_data.file_list;

for idx=1:length(file_list)

%Generate PMT file name
PMT_file=[file_list{idx} '.tdms'];

%Grab target and threshold 
target= str2num(handles.target_time.String)*1e-9;
threshold= str2num(handles.threshold_RMS.String);

%Calculate flyer delay
rise_delay=peakdelay_rising_edge_v5(PMT_file,PMT_p,target,threshold);

%Update tables with delay
handles.main_GUI_data.main_data_table{idx,4}=round(rise_delay*1e9,0);
handles.binning_table.Data{idx,2}=round(rise_delay*1e9,0);
handles.main_GUI_data.central_table.Data{idx,4}=round(rise_delay*1e9,0);


handles.delay_flag=1;
handles.rise_delay=rise_delay;
handles.rise_idx=idx;
guidata(hObject,handles);

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);

handles.delay_flag=0;
guidata(hObject,handles);

end

% --------------------------------------------------------------------
function peak_single_Callback(hObject, eventdata, handles)
% hObject    handle to peak_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;PMT_p=handles.main_GUI_data.PMT_p;
file_list=handles.main_GUI_data.file_list;

%Generate PMT file name
PMT_file=[file_list{file_idx} '.tdms'];

%Grab target and threshold 
target= str2num(handles.target_time.String)*1e-9;
threshold= str2num(handles.threshold_RMS.String);

%Calculate flyer delay
peak_delay=peakdelay_v5(PMT_file,PMT_p,target,threshold);

%Update tables with delay
handles.main_GUI_data.main_data_table{file_idx,4}=round(peak_delay*1e9,0);
handles.binning_table.Data{file_idx,2}=round(peak_delay*1e9,0);
handles.main_GUI_data.central_table.Data{file_idx,4}=round(peak_delay*1e9,0);

handles.delay_flag=2;
handles.peak_delay=peak_delay;
handles.peak_idx=file_idx;
guidata(hObject,handles);

table_single_Callback(hObject, eventdata, handles)
handles=guidata(handles.binning_GUI_bkg);

handles.delay_flag=0;
guidata(hObject,handles);

% --------------------------------------------------------------------
function peak_all_Callback(hObject, eventdata, handles)
% hObject    handle to peak_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_p=handles.main_GUI_data.PMT_p;
file_list=handles.main_GUI_data.file_list;


for idx=1:length(file_list);

%Generate PMT file name
PMT_file=[file_list{idx} '.tdms'];

%Grab target and threshold 
target= str2num(handles.target_time.String)*1e-9;
threshold= str2num(handles.threshold_RMS.String);

%Calculate flyer delay
peak_delay=peakdelay_v5(PMT_file,PMT_p,target,threshold);

%Update tables with delay
handles.main_GUI_data.main_data_table{idx,4}=round(peak_delay*1e9,0);
handles.binning_table.Data{idx,2}=round(peak_delay*1e9,0);
handles.main_GUI_data.central_table.Data{idx,4}=round(peak_delay*1e9,0);

handles.delay_flag=2;
handles.peak_delay=peak_delay;
handles.peak_idx=idx;
guidata(hObject,handles);

table_single_Callback(hObject,handles.binning_GUI_bkg,handles)
handles=guidata(handles.binning_GUI_bkg);

handles.delay_flag=0;
guidata(hObject,handles);

end



% --------------------------------------------------------------------
function table_single_Callback(hObject, eventdata, handles)
% hObject    handle to table_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;PMT_p=handles.main_GUI_data.PMT_p;
parent_p=handles.main_GUI_data.parent_p;
file_list=handles.main_GUI_data.file_list;



 %Extract data from timing table for updating
binning_data=handles.binning_table.Data;
 
 %Check for calibration lamp file
 try
    unit_conv=handles.unit_conv;
catch
cal_load_Callback(hObject, eventdata, handles);
handles=guidata(handles.binning_GUI_bkg);
unit_conv=handles.unit_conv;
 end

 
 
 switch handles.delay_flag
     case 0
 
%Obtain delay from table
%Check data table
%Check that table has a value at this point
if isempty(binning_data{file_idx,2})~=1
    delay=(cell2mat(binning_data(file_idx,2)))*1e-9;
else
    %Enter flyer delay
    prompt={'Enter impact delay (ns)'};
    dlg_title='Impact delay';
    delay=str2double(cell2mat(inputdlg(prompt,dlg_title)))*1e-9;
    
    %Update tables with delay
    handles.main_GUI_data.main_data_table{file_idx,4}=delay*1e9;
    handles.binning_table.Data{file_idx,2}=delay*1e9;
    handles.main_GUI_data.central_table.Data{file_idx,4}=delay*1e9;

end

     case 1 
         delay=handles.rise_delay;
         file_idx=handles.rise_idx;
     case 2
         delay=handles.peak_delay;
         file_idx=handles.peak_idx;
 end
 


%Generate new folder
PMT_short_p=[parent_p,'\PMT_short'];

%Generate PMT file name
PMT_file=[file_list{file_idx} '.tdms'];

%If new folder doesn't exist, create it
if exist(PMT_short_p)~=7
mkdir(PMT_short_p)
end

%translate handles
e_filter=handles.emiss_filter_flag;
c_filter=handles.cal_filter_flag;
str_dec=str2num(handles.str_dec.String);
end_dec=str2num(handles.end_dec.String);
bin_res=str2num(handles.bin_res.String);
s_chan=str2num(handles.supress_channel.String);

%shorten file
[bn_t, bn_sr]=TDMS_short_fn2_v5(PMT_file,PMT_p,PMT_short_p,delay,...
    unit_conv',e_filter,c_filter, str_dec,end_dec,bin_res,s_chan);

%Integrate spectral radiance to radiance
bn_r=sum(bn_sr,2)*1e-9;

handles.rad_plot.XData= bn_t;handles.rad_plot.YData = bn_r;

%Save data to main data cell
handles.main_GUI_data.main_data_table{file_idx,5}=bn_t;
handles.main_GUI_data.main_data_table{file_idx,6}=bn_r;
handles.main_GUI_data.main_data_table{file_idx,7}=bn_sr;

%update table with radiance
handles.main_GUI_data.central_table.Data{file_idx,7}='True';
handles.binning_table.Data{file_idx,4}='True';


for q=1:32
    handles.srad_plot(q).XData=bn_t;handles.srad_plot(q).YData=bn_sr(:,q);
end

drawnow
guidata(hObject,handles);


% --------------------------------------------------------------------
function table_all_Callback(hObject, eventdata, handles)
% hObject    handle to table_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_p=handles.main_GUI_data.PMT_p;
parent_p=handles.main_GUI_data.parent_p;
file_list=handles.main_GUI_data.file_list;

 %Extract data from timing table for updating
binning_data=handles.binning_table.Data;
 
 %Check for calibration lamp file
 try
    unit_conv=handles.unit_conv;
catch
cal_load_Callback(hObject, eventdata, handles);
handles=guidata(handles.binning_GUI_bkg);
unit_conv=handles.unit_conv;
 end

 %Loop over all files
 for idx=1:length(file_list)
 
 
%Obtain delay from table
%Check data table
%Check that table has a value at this point
if isempty(binning_data{idx,2})~=1
    delay=(cell2mat(binning_data(idx,2)))*1e-9;
else
    %Enter flyer delay
    prompt={'Enter impact delay (ns)'};
    dlg_title='Impact delay';
    delay=str2double(cell2mat(inputdlg(prompt,dlg_title)))*1e-9;

    %Update tables with delay
    handles.main_GUI_data.main_data_table{idx,4}=delay;
    handles.binning_table.Data{idx,2}=delay;
    handles.main_GUI_data.central_table.Data{idx,4}=delay*1e9;
    
end



%Generate new folder
PMT_short_p=[parent_p,'\PMT_short'];

%Generate PMT file name
PMT_file=[file_list{idx} '.tdms'];

%If new folder doesn't exist, create it
if exist(PMT_short_p)~=7
mkdir(PMT_short_p)
end

%translate handles
e_filter=handles.emiss_filter_flag;
c_filter=handles.cal_filter_flag;
str_dec=str2num(handles.str_dec.String);
end_dec=str2num(handles.end_dec.String);
bin_res=str2num(handles.bin_res.String);
s_chan=str2num(handles.supress_channel.String);

%shorten file
[bn_t, bn_sr]=TDMS_short_fn2_v5(PMT_file,PMT_p,PMT_short_p,delay,...
    unit_conv',e_filter,c_filter, str_dec,end_dec,bin_res,s_chan);

%Integrate spectral radiance to radiance
bn_r=sum(bn_sr,2)*1e-9;

handles.rad_plot.XData= bn_t;handles.rad_plot.YData = bn_r;

for q=1:32
    handles.srad_plot(q).XData=bn_t;handles.srad_plot(q).YData=bn_sr(:,q);
end

drawnow

%update main table data
handles.main_GUI_data.main_data_table{idx,5}=bn_t;
handles.main_GUI_data.main_data_table{idx,6}=bn_r;
handles.main_GUI_data.main_data_table{idx,7}=bn_sr;

%update table with radiance
handles.main_GUI_data.central_table.Data{idx,7}='True';
handles.binning_table.Data{idx,4}='True';

guidata(hObject,handles);
 end



function threshold_RMS_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_RMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_RMS as text
%        str2double(get(hObject,'String')) returns contents of threshold_RMS as a double


% --- Executes during object creation, after setting all properties.
function threshold_RMS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_RMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in binning_table.
function binning_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to binning_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
try
selected_cells=eventdata.Indices(1);

file_list=handles.main_GUI_data.file_list;

%Move to the file in the list which was selected
file_idx=selected_cells;

%Update the PMT file
PMT_file=[file_list{file_idx},'.tdms'];

%Set display
set(handles.selected_file_txt,'String',file_list{file_idx});

handles.file_idx=file_idx;
handles.PMT_file=PMT_file;
guidata(hObject,handles);
catch
end


% --- Executes on button press in cal_no_filter.
function cal_no_filter_Callback(hObject, eventdata, handles)
% hObject    handle to cal_no_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cal_no_filter.Value=1;
handles.cal03_filter.Value=0;
handles.cal06_filter.Value=0;
handles.cal13_filter.Value=0;
handles.calcube_filter.Value=0;



handles.cal_filter_flag=0;
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in cal03_filter.
function cal03_filter_Callback(hObject, eventdata, handles)
% hObject    handle to cal03_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cal_no_filter.Value=0;
handles.cal03_filter.Value=1;
handles.cal06_filter.Value=0;
handles.cal13_filter.Value=0;
handles.calcube_filter.Value=0;

handles.cal_filter_flag=1;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cal06_filter.
function cal06_filter_Callback(hObject, eventdata, handles)
% hObject    handle to cal06_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cal_no_filter.Value=0;
handles.cal03_filter.Value=0;
handles.cal06_filter.Value=1;
handles.cal13_filter.Value=0;
handles.calcube_filter.Value=0;

handles.cal_filter_flag=2;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cal13_filter.
function cal13_filter_Callback(hObject, eventdata, handles)
% hObject    handle to cal13_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cal_no_filter.Value=0;
handles.cal03_filter.Value=0;
handles.cal06_filter.Value=0;
handles.cal13_filter.Value=1;
handles.calcube_filter.Value=0;

handles.cal_filter_flag=3;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in calcube_filter.
function calcube_filter_Callback(hObject, eventdata, handles)
% hObject    handle to calcube_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cal_no_filter.Value=0;
handles.cal03_filter.Value=0;
handles.cal06_filter.Value=0;
handles.cal13_filter.Value=0;
handles.calcube_filter.Value=1;

handles.cal_filter_flag=4;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in emission_no_filter.
function emission_no_filter_Callback(hObject, eventdata, handles)
% hObject    handle to emission_no_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emission_no_filter.Value=1;
handles.emission03_filter.Value=0;
handles.emission06_filter.Value=0;
handles.emission13_filter.Value=0;
handles.emissioncube_filter.Value=0;

handles.emiss_filter_flag=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in emission03_filter.
function emission03_filter_Callback(hObject, eventdata, handles)
% hObject    handle to emission03_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emission_no_filter.Value=0;
handles.emission03_filter.Value=1;
handles.emission06_filter.Value=0;
handles.emission13_filter.Value=0;
handles.emissioncube_filter.Value=0;

handles.emiss_filter_flag=1;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in emission06_filter.
function emission06_filter_Callback(hObject, eventdata, handles)
% hObject    handle to emission06_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emission_no_filter.Value=0;
handles.emission03_filter.Value=0;
handles.emission06_filter.Value=1;
handles.emission13_filter.Value=0;
handles.emissioncube_filter.Value=0;

handles.emiss_filter_flag=2;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in emission13_filter.
function emission13_filter_Callback(hObject, eventdata, handles)
% hObject    handle to emission13_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emission_no_filter.Value=0;
handles.emission03_filter.Value=0;
handles.emission06_filter.Value=0;
handles.emission13_filter.Value=1;
handles.emissioncube_filter.Value=0;

handles.emiss_filter_flag=3;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in emissioncube_filter.
function emissioncube_filter_Callback(hObject, eventdata, handles)
% hObject    handle to emissioncube_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emission_no_filter.Value=0;
handles.emission03_filter.Value=0;
handles.emission06_filter.Value=0;
handles.emission13_filter.Value=0;
handles.emissioncube_filter.Value=1;

handles.emiss_filter_flag=4;
% Update handles structure
guidata(hObject, handles);



function bin_res_Callback(hObject, eventdata, handles)
% hObject    handle to bin_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function bin_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bin_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function str_dec_Callback(hObject, eventdata, handles)
% hObject    handle to str_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function str_dec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to str_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_dec_Callback(hObject, eventdata, handles)
% hObject    handle to end_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function end_dec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function cal_load_Callback(hObject, eventdata, handles)
% hObject    handle to cal_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Assign center wavelength for fiber
wavelength=[444.88,450.88,457.38,464.38,471.88,480,488.88,498.5,509,520.5,...
    530.88,540,549.75,560.13,571.38,583.63,596.88,611.25,623.75,634,644.88,...
    656.5,668.88,682,696,711,727.13,744.5,763.25,783.5,805.38,829.25]';

%Import data
raw_file=convertTDMS_cal(0);

%Extract data
A=extractfield(raw_file.Data.MeasuredData,'Data');

%Shape data
B=reshape(A,length(A)/33,33);

%Isolate time and flip voltage
time=B(:,5)+B(1,5);B(:,5)=[];
fiber_array=B*-1;

%Move from volts to current at 1 Mohm terminator
PMT_current=fiber_array*1e-6;

%Average the flat current
avg_current=mean(PMT_current,1)';

handles.radcal_figure=figure('Name','Calibration with Radiance Lamp',...
    'Tag','radiance calibration figure','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);

%plot raw radiance in current
handles.radcal_raw=subplot(2,1,1);plot(time,PMT_current);
set(get(handles.radcal_raw, 'yLabel'), 'String', 'current (A)');

%Divide by radiation lamp values
lamp=[80037000,90076000,98599000,108250000,120900000,134170000,...
    150670000,169850000,193110000,218410000,245070000,270140000,...
    300020000,333850000,370060000,407710000,451600000,487930000,...
    520190000,548530000,585790000,613660000,624730000,648610000,...
    683880000,710760000,740110000,776370000,812140000,847370000,...
    883530000,923480000]';


%Determine Calibration vector
unit_conversion=lamp./avg_current.*0.1;

assignin('base','unit_conversion',unit_conversion);

%Plot calibration vector
handles.radcal_vector=subplot(2,1,2);plot(wavelength,unit_conversion,'Marker','o');
set(get(handles.radcal_vector, 'yLabel'), 'String', 'W/sr \cdot m^3 \cdot A');
set(get(handles.radcal_vector, 'xLabel'), 'String', 'wavelength (nm)');

%Update table with calibration
handles.main_GUI_data.main_data_table{1,15}=unit_conversion;

handles.unit_conv=unit_conversion;

guidata(hObject,handles);



function supress_channel_Callback(hObject, eventdata, handles)
% hObject    handle to supress_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of supress_channel as text
%        str2double(get(hObject,'String')) returns contents of supress_channel as a double


% --- Executes during object creation, after setting all properties.
function supress_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to supress_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in left_refine.
function left_refine_Callback(hObject, eventdata, handles)
% hObject    handle to left_refine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;
file_list=handles.main_GUI_data.file_list;

%Move to the previous file in the list if there is one
if file_idx-1 > 0;
    file_idx=file_idx-1;
end

%Set display
set(handles.selected_file_txt,'String',file_list{file_idx});

handles.file_idx=file_idx;
guidata(hObject,handles);


table_single_Callback(hObject, eventdata, handles)


% --- Executes on button press in right_refine.
function right_refine_Callback(hObject, eventdata, handles)
% hObject    handle to right_refine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;max_length=handles.main_GUI_data.max_length;
file_list=handles.main_GUI_data.file_list;

%Move to the next file in the list if there is one
if file_idx+1 <= max_length;
    file_idx=file_idx+1;
end

%Set display
set(handles.selected_file_txt,'String',file_list{file_idx});

handles.file_idx=file_idx;
guidata(hObject,handles);

table_single_Callback(hObject, eventdata, handles)


% --- Executes when entered data in editable cell(s) in binning_table.
function binning_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to binning_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%Identify changed value
delay=eventdata.NewData;
cell=eventdata.Indices(1);

handles.binning_table.Data{cell,2}=delay;
handles.main_GUI_data.main_data_table{cell,4}=delay;

guidata(hObject,handles);


% --- Executes when user attempts to close binning_GUI_bkg.
function binning_GUI_bkg_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to binning_GUI_bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.decision]=GUI_pass;

switch handles.decision
    case 1
        main_GUI=findobj('Tag','Emission_spectrometer_GUI');
        guidata(main_GUI,handles.main_GUI_data);
        delete(handles.binning_GUI_bkg);
    case 2
        delete(handles.binning_GUI_bkg);
    case 3
end


% --------------------------------------------------------------------
function pass_close_Callback(hObject, eventdata, handles)
% hObject    handle to pass_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.decision=1;
guidata(hObject,handles);

main_GUI=findobj('Tag','Emission_spectrometer_GUI');
guidata(main_GUI,handles.main_GUI_data);
delete(handles.binning_GUI_bkg);
