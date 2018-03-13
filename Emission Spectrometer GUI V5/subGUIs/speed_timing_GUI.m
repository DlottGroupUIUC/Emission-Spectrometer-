function varargout = speed_timing_GUI(varargin)
% SPEED_TIMING_GUI MATLAB code for speed_timing_GUI.fig
%      SPEED_TIMING_GUI, by itself, creates a new SPEED_TIMING_GUI or raises the existing
%      singleton*.
%
%      H = SPEED_TIMING_GUI returns the handle to a new SPEED_TIMING_GUI or the handle to
%      the existing singleton*.
%
%      SPEED_TIMING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPEED_TIMING_GUI.M with the given input arguments.
%
%      SPEED_TIMING_GUI('Property','Value',...) creates a new SPEED_TIMING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before speed_timing_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to speed_timing_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help speed_timing_GUI

% Last Modified by GUIDE v2.5 13-Mar-2018 13:53:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @speed_timing_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @speed_timing_GUI_OutputFcn, ...
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


% --- Executes just before speed_timing_GUI is made visible.
function speed_timing_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to speed_timing_GUI (see VARARGIN)
%load data from main GUI
h = findobj('Tag','Emission_spectrometer_GUI');
 
  % if exists (not empty)
 if ~isempty(h)
    % get handles and other user-defined data associated to Gui1
    main_GUI_data = guidata(h);

 end
 
%Set file index
file_idx=1;

%Set table data
handles.S_T_table.Data=main_GUI_data.central_table.Data(:,[1:4,6]);

%Set current plot
handles.selected_file_txt.String=main_GUI_data.file_list{file_idx};

%%%Initialize Plots%%%
%Generate data
time=linspace(-500,1500,1e4);
vel=zeros(1,1e4);
volts=zeros(1,1e4);
vel_img=zeros(1e4,1e4);

handles.vel_sub=subplot('Position',[0.28 0.68 0.70 0.30]);handles.vel_plot=plot(time,vel,'b',time,vel,'r',time,vel,'k');
handles.int_sub=subplot('Position',[0.28 0.28 0.70 0.10]);handles.int_plot=plot(time,volts,'b',time,volts,'r',time,volts,'k');
handles.spec_sub=subplot('Position',[0.28 0.38 0.70 0.30]);handles.spec_plot=imagesc(time,vel,vel_img);
handles.spec_LO_sub=subplot('Position',[0.28,0.08,0.70,0.15]);handles.spec_LO_plot=plot(vel,volts,'k',vel,volts,'r',vel,volts,'b');


handles.spec_sub.XTickLabel=[];
handles.vel_sub.XTickLabel=[];
ylabel(handles.vel_sub,'velocity (km/s)');
ylabel(handles.int_sub,'volts (V)');
xlabel(handles.int_sub,'time (s)');
ylabel(handles.spec_sub,'velocity (km/s)')
linkaxes([handles.vel_sub handles.int_sub, handles.spec_sub],'x');
linkaxes([handles.vel_sub handles.spec_sub],'y');
set(handles.spec_sub,'ydir','normal')

ylabel(handles.spec_LO_sub,'intensity (au)');
xlabel(handles.spec_LO_sub,'velocity (km/s)');
handles.spec_LO_sub.XLimMode='auto';
handles.spec_LO_sub.YLimMode='auto';
handles.vel_sub.XLimMode='auto';
handles.vel_sub.YLimMode='auto';

%Initialize PDV vertical cutoff
handles.PDV_v_cutoff.String='0.1';
handles.X_min_tag.String='-500';
handles.X_max_tag.String='1500';
handles.Y_max_tag.String='6';

% Update main GUI
handles.output = hObject;
handles.main_GUI_data=main_GUI_data;
handles.file_idx=file_idx;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes speed_timing_GUI wait for user response (see UIRESUME)
% uiwait(handles.speed_timing_GUI_fig);


% --- Outputs from this function are returned to the command line.
function varargout = speed_timing_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function S_T_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to S_T_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%initialize table
set(hObject,'Data',cell(1,5));
set(hObject, 'ColumnName', {'shot #','speed (km/s)','error (km/s)','impact time (ns)','PDV' });
set(hObject, 'ColumnFormat',{'numeric','numeric','numeric','numeric','logical'});


% --------------------------------------------------------------------
function speed_header_Callback(hObject, eventdata, handles)
% hObject    handle to speed_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function impact_header_Callback(hObject, eventdata, handles)
% hObject    handle to impact_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function impact_single_Callback(hObject, eventdata, handles)
% hObject    handle to impact_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;PDV_p=handles.main_GUI_data.PDV_p;
file_list=handles.main_GUI_data.file_list;

%Determine if a file exists here
if isempty(handles.S_T_table.Data{file_idx,5})~=1
    
    %Update current shot
    set(handles.selected_file_txt,'String',file_list{file_idx});
    
    %Run full time analysis for channel 1
    [t,s,time,intensity]=fulltime_fn_v5(file_list{file_idx},PDV_p,1);
    %Run full time analysis for channel 2
    [t2,s2,time2,intensity2]=fulltime_fn_v5(file_list{file_idx},PDV_p,2);
    %Run full time analysis for channel 4
    [t4,s4,time4,intensity4]=fulltime_fn_v5(file_list{file_idx},PDV_p,4);
    
    %combine channels
    combined=sortrows([t,s;t2,s2;t4,s4]);
    combined_time=combined(:,1);
    combined_speed=combined(:,2);
    
    %Update plots before asking for cursors
    handles.vel_plot(1).XData=t;handles.vel_plot(1).YData=s;
    handles.vel_plot(2).XData=t2;handles.vel_plot(2).YData=s2;
    handles.vel_plot(3).XData=t4;handles.vel_plot(3).YData=s4;
    
    handles.int_plot(1).XData=time;handles.int_plot(1).YData=intensity;
    handles.int_plot(2).XData=time2;handles.int_plot(2).YData=intensity2;
    handles.int_plot(3).XData=time4;handles.int_plot(3).YData=intensity4;
    
    handles.spec_plot.XData=t;
    handles.spec_plot.YData=zeros(1,length(t));
    handles.spec_plot.CData=zeros(length(t));
    
    %Select impact time
    [impact_time ~]=ginput_crosshair_v4(1);
    
    %update plots
    handles.spec_LO_plot(1).XData=zeros(1,1e3);handles.spec_LO_plot(1).YData=zeros(1,1e3);
    handles.spec_LO_plot(2).XData=zeros(1,1e3);handles.spec_LO_plot(2).YData=zeros(1,1e3);
    handles.spec_LO_plot(3).XData=zeros(1,1e3);handles.spec_LO_plot(3).YData=zeros(1,1e3);
    
    
    %Update tables
    handles.main_GUI_data.main_data_table{file_idx,4}=round(impact_time);
    
    handles.S_T_table.Data{file_idx,4}=round(impact_time);
    
    handles.main_GUI_data.central_table.Data{file_idx,4}=round(impact_time);
    
end

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function impact_all_Callback(hObject, eventdata, handles)
% hObject    handle to impact_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PDV_p=handles.main_GUI_data.PDV_p;
file_list=handles.main_GUI_data.file_list;

%loop files
for idx=1:length(file_list)
    
    %Determine if a file exists here
    if isempty(handles.S_T_table.Data{idx,5})~=1
        
        %Update current shot
        set(handles.selected_file_txt,'String',file_list{idx});
        
        %Run full time analysis for channel 1
        [t,s,time,intensity]=fulltime_fn_v5(file_list{idx},PDV_p,1);
        %Run full time analysis for channel 2
        [t2,s2,time2,intensity2]=fulltime_fn_v5(file_list{idx},PDV_p,2);
        %Run full time analysis for channel 4
        [t4,s4,time4,intensity4]=fulltime_fn_v5(file_list{idx},PDV_p,4);
        
        %combine channels
        combined=sortrows([t,s;t2,s2;t4,s4]);
        combined_time=combined(:,1);
        combined_speed=combined(:,2);
        
        %Update plots before asking for cursors
        handles.vel_plot(1).XData=t;handles.vel_plot(1).YData=s;
        handles.vel_plot(2).XData=t2;handles.vel_plot(2).YData=s2;
        handles.vel_plot(3).XData=t4;handles.vel_plot(3).YData=s4;
        
        handles.int_plot(1).XData=time;handles.int_plot(1).YData=intensity;
        handles.int_plot(2).XData=time2;handles.int_plot(2).YData=intensity2;
        handles.int_plot(3).XData=time4;handles.int_plot(3).YData=intensity4;
        
        handles.spec_plot.XData=t;
        handles.spec_plot.YData=zeros(1,length(t));
        handles.spec_plot.CData=zeros(length(t));
        
        %Select impact time
        [impact_time ~]=ginput_crosshair_v4(1);
        
        %update plots
        handles.spec_LO_plot(1).XData=zeros(1,1e3);handles.spec_LO_plot(1).YData=zeros(1,1e3);
        handles.spec_LO_plot(2).XData=zeros(1,1e3);handles.spec_LO_plot(2).YData=zeros(1,1e3);
        handles.spec_LO_plot(3).XData=zeros(1,1e3);handles.spec_LO_plot(3).YData=zeros(1,1e3);
        
        
        %Update tables
        handles.main_GUI_data.main_data_table{idx,4}=round(impact_time);
        
        handles.S_T_table.Data{idx,4}=round(impact_time);
        
        handles.main_GUI_data.central_table.Data{idx,4}=round(impact_time);
        
    end
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function speed_single_Callback(hObject, eventdata, handles)
% hObject    handle to speed_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;PDV_p=handles.main_GUI_data.PDV_p;
file_list=handles.main_GUI_data.file_list;
PDV_v_cutoff=str2num(handles.PDV_v_cutoff.String);

%Determine if a file exists here
if isempty(handles.S_T_table.Data{file_idx,5})~=1
    
    %Update current shot
    set(handles.selected_file_txt,'String',file_list{file_idx});
    
    %Run STFT for channel 1
    [t,s,time,intensity,SG_v1,SG_t1,SG_1]=STFT_fn_v5(file_list{file_idx},PDV_p,1,PDV_v_cutoff);
    %Run STFT time for channel 2
    [t2,s2,time2,intensity2,SG_v2,~,SG_2]=STFT_fn_v5(file_list{file_idx},PDV_p,2,PDV_v_cutoff);
    %Run STFT time for channel 4
    [t4,s4,time4,intensity4,SG_v4,~,SG_4]=STFT_fn_v5(file_list{file_idx},PDV_p,4,PDV_v_cutoff);
    
    %combine line outs
    combined=sortrows([t,s;t2,s2;t4,s4]);
    combined_time=combined(:,1);
    combined_speed=combined(:,2);
    
    %combine spectrograms
    comb_SGv=[SG_v1;SG_v2];
    comb_SG=[SG_1;SG_2];
    sort_SG=sortrows([comb_SGv,comb_SG]);
    sorted_SGv=sort_SG(:,1);
    sorted_SGSG=sort_SG(:,[2:end]);
    
    %Update plots before asking for cursors
    handles.vel_plot(1).XData=t;handles.vel_plot(1).YData=s;
    handles.vel_plot(2).XData=t2;handles.vel_plot(2).YData=s2;
    handles.vel_plot(3).XData=t4;handles.vel_plot(3).YData=s4;
    
    handles.int_plot(1).XData=time;handles.int_plot(1).YData=intensity;
    handles.int_plot(2).XData=time2;handles.int_plot(2).YData=intensity2;
    handles.int_plot(3).XData=time4;handles.int_plot(3).YData=intensity4;
    
    handles.spec_plot.XData=SG_t1;
    handles.spec_plot.YData=sorted_SGv;
    handles.spec_plot.CData=sorted_SGSG;
    
    %Grab points for LTFT
    [xm,ym] = ginput_crosshair_v4(1); %xmouse, ymouse
    [~, yidx] = min(abs(combined_speed-ym)); %closest index y
    [~, xidx] = min(abs(combined_time-xm)); %closest index x
    theTime1=combined_time(xidx); %extract x
    
    [xm,ym] = ginput_crosshair_v4(1); %xmouse, ymouse
    [~, yidx] = min(abs(combined_speed-ym)); %closest index y
    [~, xidx] = min(abs(combined_time-xm)); %closest index x
    theSpeed = combined_speed(yidx); %extract y
    theTime2=combined_time(xidx); %extract x
    
    [velocity,unc,SG_line_v1,SG_line_FFT1]=LTFT_error_fn_v5(time,intensity,theTime1,theTime2,theSpeed);
    [velocity2,unc2,SG_line_v2,SG_line_FFT2]=LTFT_error_fn_v5(time2,intensity2,theTime1,theTime2,theSpeed);
    [velocity4,unc4,SG_line_v4,SG_line_FFT4]=LTFT_error_fn_v5(time4,intensity4,theTime1,theTime2,theSpeed);
    
    %update plots
    handles.spec_LO_plot(1).XData=SG_line_v1;handles.spec_LO_plot(1).YData=SG_line_FFT1;
    handles.spec_LO_plot(2).XData=SG_line_v2;handles.spec_LO_plot(2).YData=SG_line_FFT2;
    handles.spec_LO_plot(3).XData=SG_line_v4;handles.spec_LO_plot(3).YData=SG_line_FFT4;
    
    %Update tables
    handles.main_GUI_data.main_data_table{file_idx,2}=round(mean([velocity velocity2 velocity4])*100)/100;
    handles.main_GUI_data.main_data_table{file_idx,3}=mean([unc unc2 unc4]);
    
    handles.S_T_table.Data{file_idx,2}=round(mean([velocity velocity2 velocity4])*100)/100;
    handles.S_T_table.Data{file_idx,3}=mean([unc unc2 unc4]);
    
    handles.main_GUI_data.central_table.Data{file_idx,2}=round(mean([velocity velocity2 velocity4])*100)/100;
    handles.main_GUI_data.central_table.Data{file_idx,3}=mean([unc unc2 unc4]);
    
end

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function speed_all_Callback(hObject, eventdata, handles)
% hObject    handle to speed_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PDV_p=handles.main_GUI_data.PDV_p;
file_list=handles.main_GUI_data.file_list;
PDV_v_cutoff=str2num(handles.PDV_v_cutoff.String);

for idx=1:length(file_list);
    
    %Determine if a file exists here
    if isempty(handles.S_T_table.Data{file_idx,5})~=1
        
        %Update current shot
        set(handles.selected_file_txt,'String',file_list{idx});
        
        %Run STFT for channel 1
        [t,s,time,intensity,SG_v1,SG_t1,SG_1]=STFT_fn_v5(file_list{idx},PDV_p,1,PDV_v_cutoff);
        %Run STFT time for channel 2
        [t2,s2,time2,intensity2,SG_v2,~,SG_2]=STFT_fn_v5(file_list{idx},PDV_p,2,PDV_v_cutoff);
        %Run STFT time for channel 4
        [t4,s4,time4,intensity4,SG_v4,~,SG_4]=STFT_fn_v5(file_list{idx},PDV_p,4,PDV_v_cutoff);
        
        %combine line outs
        combined=sortrows([t,s;t2,s2;t4,s4]);
        combined_time=combined(:,1);
        combined_speed=combined(:,2);
        
        %combine spectrograms
        comb_SGv=[SG_v1;SG_v2];
        comb_SG=[SG_1;SG_2];
        sort_SG=sortrows([comb_SGv,comb_SG]);
        sorted_SGv=sort_SG(:,1);
        sorted_SGSG=sort_SG(:,[2:end]);
        
        %Update plots before asking for cursors
        handles.vel_plot(1).XData=t;handles.vel_plot(1).YData=s;
        handles.vel_plot(2).XData=t2;handles.vel_plot(2).YData=s2;
        handles.vel_plot(3).XData=t4;handles.vel_plot(3).YData=s4;
        
        handles.int_plot(1).XData=time;handles.int_plot(1).YData=intensity;
        handles.int_plot(2).XData=time2;handles.int_plot(2).YData=intensity2;
        handles.int_plot(3).XData=time4;handles.int_plot(3).YData=intensity4;
        
        handles.spec_plot.XData=SG_t1;
        handles.spec_plot.YData=sorted_SGv;
        handles.spec_plot.CData=sorted_SGSG;
        
        %Grab points for LTFT
        [xm,ym] = ginput_crosshair_v4(1); %xmouse, ymouse
        [~, yidx] = min(abs(combined_speed-ym)); %closest index y
        [~, xidx] = min(abs(combined_time-xm)); %closest index x
        theTime1=combined_time(xidx); %extract x
        
        [xm,ym] = ginput_crosshair_v4(1); %xmouse, ymouse
        [~, yidx] = min(abs(combined_speed-ym)); %closest index y
        [~, xidx] = min(abs(combined_time-xm)); %closest index x
        theSpeed = combined_speed(yidx); %extract y
        theTime2=combined_time(xidx); %extract x
        
        [velocity,unc,SG_line_v1,SG_line_FFT1]=LTFT_error_fn_v5(time,intensity,theTime1,theTime2,theSpeed);
        [velocity2,unc2,SG_line_v2,SG_line_FFT2]=LTFT_error_fn_v5(time2,intensity2,theTime1,theTime2,theSpeed);
        [velocity4,unc4,SG_line_v4,SG_line_FFT4]=LTFT_error_fn_v5(time4,intensity4,theTime1,theTime2,theSpeed);
        
        %update plots
        handles.spec_LO_plot(1).XData=SG_line_v1;handles.spec_LO_plot(1).YData=SG_line_FFT1;
        handles.spec_LO_plot(2).XData=SG_line_v2;handles.spec_LO_plot(2).YData=SG_line_FFT2;
        handles.spec_LO_plot(3).XData=SG_line_v4;handles.spec_LO_plot(3).YData=SG_line_FFT4;
        
        %Update tables
        handles.main_GUI_data.main_data_table{idx,2}=round(mean([velocity velocity2 velocity4])*100)/100;
        handles.main_GUI_data.main_data_table{idx,3}=mean([unc unc2 unc4]);
        
        handles.S_T_table.Data{idx,2}=round(mean([velocity velocity2 velocity4])*100)/100;
        handles.S_T_table.Data{idx,3}=mean([unc unc2 unc4]);
        
        handles.main_GUI_data.central_table.Data{idx,2}=round(mean([velocity velocity2 velocity4])*100)/100;
        handles.main_GUI_data.central_table.Data{idx,3}=mean([unc unc2 unc4]);
        
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function selected_file_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_file_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in spectogram_flag.
function spectogram_flag_Callback(hObject, eventdata, handles)
% hObject    handle to spectogram_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function pass_close_Callback(hObject, eventdata, handles)
% hObject    handle to pass_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.decision=1;
guidata(hObject,handles);

main_GUI=findobj('Tag','Emission_spectrometer_GUI');
guidata(main_GUI,handles.main_GUI_data);
delete(handles.speed_timing_GUI_fig);


% --- Executes when selected cell(s) is changed in S_T_table.
function S_T_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to S_T_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
try
selected_cells=eventdata.Indices(1);

file_list=handles.main_GUI_data.file_list;

%Move to the file in the list which was selected
file_idx=selected_cells;

%Update the PMT file
PDV_file=[file_list{file_idx},'.tdms'];

%Set display
set(handles.selected_file_txt,'String',file_list{file_idx});

handles.file_idx=file_idx;
handles.PDV_file=PDV_file;
guidata(hObject,handles);
catch
end


% --- Executes when user attempts to close speed_timing_GUI_fig.
function speed_timing_GUI_fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to speed_timing_GUI_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.decision]=GUI_pass;

switch handles.decision
    case 1
        main_GUI=findobj('Tag','Emission_spectrometer_GUI');
        guidata(main_GUI,handles.main_GUI_data);
        delete(handles.speed_timing_GUI_fig);
    case 2
        delete(handles.speed_timing_GUI_fig);
    case 3
end



function PDV_v_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to PDV_v_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PDV_v_cutoff as text
%        str2double(get(hObject,'String')) returns contents of PDV_v_cutoff as a double


% --- Executes during object creation, after setting all properties.
function PDV_v_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PDV_v_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function X_max_tag_Callback(hObject, eventdata, handles)
% hObject    handle to X_max_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.vel_sub.XLim(1)=str2num(handles.X_min_tag.String);
handles.vel_sub.XLim(2)=str2num(handles.X_max_tag.String);


% --- Executes during object creation, after setting all properties.
function X_max_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_max_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function X_min_tag_Callback(hObject, eventdata, handles)
% hObject    handle to X_min_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.vel_sub.XLim(1)=str2num(handles.X_min_tag.String);
handles.vel_sub.XLim(2)=str2num(handles.X_max_tag.String);


% --- Executes during object creation, after setting all properties.
function X_min_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_min_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_max_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Y_max_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.vel_sub.YLim(2)=str2num(handles.Y_max_tag.String);


% --- Executes during object creation, after setting all properties.
function Y_max_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_max_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
