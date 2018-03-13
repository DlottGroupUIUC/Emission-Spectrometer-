function varargout = emission_processing_GUI(varargin)
% EMISSION_PROCESSING_GUI MATLAB code for emission_processing_GUI.fig
%      EMISSION_PROCESSING_GUI, by itself, creates a new EMISSION_PROCESSING_GUI or raises the existing
%      singleton*.
%
%      H = EMISSION_PROCESSING_GUI returns the handle to a new EMISSION_PROCESSING_GUI or the handle to
%      the existing singleton*.
%
%      EMISSION_PROCESSING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMISSION_PROCESSING_GUI.M with the given input arguments.
%
%      EMISSION_PROCESSING_GUI('Property','Value',...) creates a new EMISSION_PROCESSING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before emission_processing_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to emission_processing_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help emission_processing_GUI

% Last Modified by GUIDE v2.5 12-Mar-2018 10:22:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @emission_processing_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @emission_processing_GUI_OutputFcn, ...
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


% --- Executes just before emission_processing_GUI is made visible.
function emission_processing_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to emission_processing_GUI (see VARARGIN)

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

%Generate file name
PMT_short_file=[file_list{file_idx} '_short.txt.'];

%Set table data
handles.data_table.Data=main_GUI_data.central_table.Data(:,[1,7,8,9]);

%Assign center wavelength for fiber and calculate inverse wavelengths
wavelength=[444.88,450.88,457.38,464.38,471.88,480,488.88,498.5,509,520.5,...
    530.88,540,549.75,560.13,571.38,583.63,596.88,611.25,623.75,634,644.88,...
    656.5,668.88,682,696,711,727.13,744.5,763.25,783.5,805.38,829.25]';
inv_wavelength=1./(wavelength*1e-9);

%Set current plot
handles.selected_file_txt.String=main_GUI_data.file_list{file_idx};

%%%Initialize Plots%%%
%Generate data
time=logspace(-9,-5,1e3);
rad=zeros(1,1e3);
temp=zeros(1,1e3);

axes(handles.emiss_axis);
handles.temp_sub=subplot('Position',[0.25 0.58 0.70 0.40]);handles.temp_plot=semilogx(time,temp,'ok',time,temp,'dr');
handles.rad_sub=subplot('Position',[0.25 0.18 0.70 0.40]);handles.rad_plot=semilogx(time,rad,'b');

set(handles.temp_sub,'XTickLabel',[]);
ylabel(handles.rad_sub,'radiance (W/sr \cdot m^2)');
ylabel(handles.temp_sub,'temperature (K)');
xlabel(handles.rad_sub,'time (s)');
xlim(handles.rad_sub,[1e-9 1e-5]);
ylim(handles.temp_sub,[0 7000]);
linkaxes([handles.rad_sub handles.temp_sub],'x');
legend(handles.temp_sub,'graybody','Z','Location','northeast');


%initialize close function
handles.decision=3;

% Update main GUI
handles.output = hObject;
handles.main_GUI_data=main_GUI_data;
handles.file_idx=file_idx;
handles.PMT_short_file=PMT_short_file;
handles.wavelength=wavelength;
handles.inv_wavelength=inv_wavelength;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes emission_processing_GUI wait for user response (see UIRESUME)
% uiwait(handles.emission_processing_GUI_bkg);


% --- Outputs from this function are returned to the command line.
function varargout = emission_processing_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function data_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function WL_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WL_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Data',{'True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True','True'});
set(hObject, 'ColumnName',{'444.88 (1)','450.88 (2)','457.38 (3)','464.38 (4)','471.88 (5)','480 (6)','488.88 (7)','498.5 (8)','509 (9)','520.5 (10)',...
    '530.88 (11)','540 (12)','549.75 (13)','560.13 (14)','571.38 (15)','583.63 (16)','596.88 (17)','611.25 (18)','623.75 (19)','634 (20)','644.88 (21)',...
    '656.5 (22)','668.88 (23)','682 (24)','696 (25)','711 (26)','727.13 (27)','744.5 (28)','763.25 (29)','783.5 (30)','805.38 (31)','829.25 (32)'});
set(hObject, 'ColumnFormat',{'logical','logical','logical','logical',...
    'logical','logical','logical','logical','logical','logical','logical',...
    'logical','logical','logical','logical','logical','logical','logical',...
    'logical','logical','logical','logical','logical','logical','logical',...
    'logical','logical','logical','logical','logical','logical','logical'});


% --- Executes when selected cell(s) is changed in WL_table.
function WL_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to WL_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
try
    selected_cells=eventdata.Indices(2);
    
    %Toggle active state of the selected wavelength
    if strcmp(eventdata.Source.Data{selected_cells},'True')==1;
        
        eventdata.Source.Data{selected_cells}='False';
        
    elseif strcmp(eventdata.Source.Data{selected_cells},'False')==1;
        
        eventdata.Source.Data{selected_cells}='True';
        
    else
    end
catch
end


% --- Executes when selected cell(s) is changed in data_table.
function data_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to data_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
file_list=handles.main_GUI_data.file_list;

try
    selected_cells=eventdata.Indices(1);
    
    %Move to the file in the list which was selected
    file_idx=selected_cells;
    
    %Generate file name
    PMT_short_file=[file_list{file_idx} '_short.txt.'];
    
    %Set display
    set(handles.selected_file_txt,'String',file_list{file_idx});
    
    handles.file_idx=file_idx;
    handles.PMT_short_file=PMT_short_file;
    guidata(hObject,handles);
catch
end


% --------------------------------------------------------------------
function Z_header_Callback(hObject, eventdata, handles)
% hObject    handle to Z_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gray_header_Callback(hObject, eventdata, handles)
% hObject    handle to gray_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gray_single_Callback(hObject, eventdata, handles)
% hObject    handle to gray_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;
file_list=handles.main_GUI_data.file_list;wavelength=handles.wavelength;


%Check data table
table_data=get(handles.data_table,'Data');

%check for Z data. If not present, skip
if table_data{file_idx,3} == 'True'
    
    %Update current shot
    set(handles.selected_file_txt,'String',file_list{file_idx});
    
    %import necessary data
    time=handles.main_GUI_data.main_data_table{file_idx,5};
    radiance=handles.main_GUI_data.main_data_table{file_idx,6};
    spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7};
    T=handles.main_GUI_data.main_data_table{file_idx,8};
    
    
    
    %Turn off warning about negative data
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    
    %Get list of wavelengths to exclude
    wl_list=get(handles.WL_table,'Data');
    wl_exclude=find(not(cellfun('isempty', strfind(wl_list,'False'))));
    
    %Generate spectral radiance fits
    for spec_idx=1:size(spectral_radiance,1)
        if isnan(T(spec_idx))~=1
            [model(spec_idx,:),spec_T(spec_idx),spec_T_unc(spec_idx),spec_emiss(spec_idx),...
                spec_emiss_unc(spec_idx)]=spec_rad_temp_fit2_v5(wavelength,spectral_radiance(spec_idx,:),...
                T(spec_idx),wl_exclude);
        else
            spec_T(spec_idx)=NaN;spec_emiss(spec_idx)=NaN;
            spec_T_unc(spec_idx)=NaN;spec_emiss_unc(spec_idx)=NaN;
        end
    end
    
    %Eliminate fits where the uncertainty is greater than 10%
    spec_emiss(spec_T_unc>.1*spec_T)=NaN;
    spec_T(spec_T_unc>.1*spec_T)=NaN;
    
    
    %Update tables
    handles.data_table.Data{file_idx,4}='True';
    handles.main_GUI_data.central_table.Data{file_idx,9}='True';
    
    %update plots
    handles.rad_plot.XData= time;handles.rad_plot.YData = radiance;
    handles.temp_plot(1).XData= time;handles.temp_plot(1).YData = spec_T;
    handles.temp_plot(2).XData= time;handles.temp_plot(2).YData = zeros(1,length(time));
    
    %Update main data cell
    handles.main_GUI_data.main_data_table{file_idx,10}=model;
    handles.main_GUI_data.main_data_table{file_idx,11}=spec_T';
    handles.main_GUI_data.main_data_table{file_idx,12}=spec_T_unc';
    handles.main_GUI_data.main_data_table{file_idx,13}=spec_emiss';
    handles.main_GUI_data.main_data_table{file_idx,14}=spec_emiss_unc';
    
    
end


guidata(hObject,handles);


% --------------------------------------------------------------------
function gray_all_Callback(hObject, eventdata, handles)
% hObject    handle to gray_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_list=handles.main_GUI_data.file_list;wavelength=handles.wavelength;

%Check data table
table_data=get(handles.data_table,'Data');

%Run loop
for idx=1:length(file_list);
    
    
    
    %Check for Z data. If not present, skip this file
    if table_data{idx,3} == 'True';
        
        %Update current shot
        set(handles.selected_file_txt,'String',file_list{idx});
        
        
        
        
        %Import data
        time=handles.main_GUI_data.main_data_table{idx,5};
        radiance=handles.main_GUI_data.main_data_table{idx,6};
        spectral_radiance=handles.main_GUI_data.main_data_table{idx,7};
        T=handles.main_GUI_data.main_data_table{idx,8};
        
        %Turn off warning about negative data
        warning('off','MATLAB:Axes:NegativeDataInLogAxis');
        
        %Get list of wavelengths to exclude
        wl_list=get(handles.WL_table,'Data');
        wl_exclude=find(not(cellfun('isempty', strfind(wl_list,'False'))));
        
        %Generate spectral radiance fits
        parfor spec_idx=1:size(spectral_radiance,1)
            if isnan(T(spec_idx))~=1
                [model(spec_idx,:),spec_T(spec_idx),spec_T_unc(spec_idx),spec_emiss(spec_idx),...
                    spec_emiss_unc(spec_idx)]=spec_rad_temp_fit2_v5(wavelength,spectral_radiance(spec_idx,:),...
                    T(spec_idx),wl_exclude);
            else
                spec_T(spec_idx)=NaN;spec_emiss(spec_idx)=NaN;
                spec_T_unc(spec_idx)=NaN;spec_emiss_unc(spec_idx)=NaN;
            end
        end
        
        %Eliminate fits where the uncertainty is greater than 10%
        spec_emiss(spec_T_unc>.1*spec_T)=NaN;
        spec_T(spec_T_unc>.1*spec_T)=NaN;
        
        
        
        %Update tables
        handles.data_table.Data{idx,4}='True';
        handles.main_GUI_data.central_table.Data{idx,9}='True';
        
        %update plots
        handles.rad_plot.XData= time;handles.rad_plot.YData = radiance;
        handles.temp_plot(1).XData= time;handles.temp_plot(1).YData = spec_T;
        handles.temp_plot(2).XData= time;handles.temp_plot(2).YData = zeros(1,length(time));
        
        %Update main data cell
        handles.main_GUI_data.main_data_table{idx,10}=model;
        handles.main_GUI_data.main_data_table{idx,11}=spec_T;
        handles.main_GUI_data.main_data_table{idx,12}=spec_T_unc;
        handles.main_GUI_data.main_data_table{idx,13}=spec_emiss;
        handles.main_GUI_data.main_data_table{idx,14}=spec_emiss_unc;
        
        %clear variables to prevent errors
        clear time radiance spec_T spec_T_unc spec_emiss spec_emiss_unc
        drawnow
    end
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function Z_single_Callback(hObject, eventdata, handles)
% hObject    handle to Z_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;
file_list=handles.main_GUI_data.file_list;wavelength=handles.wavelength;


%Check data table
table_data=get(handles.data_table,'Data');

%Check for radiance file. If not present, skip
if table_data{file_idx,2}=='True'
    
    
    %Update current shot
    handles.selected_file_txt.String=file_list{file_idx};
    
       
    %import necessary data
    time=handles.main_GUI_data.main_data_table{file_idx,5};
    radiance=handles.main_GUI_data.main_data_table{file_idx,6};
    spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7};
    
    
    %Turn off warning about negative data
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    
    %Get list of wavelengths to exclude
    wl_list=get(handles.WL_table,'Data');
    wl_exclude=find(not(cellfun('isempty', strfind(wl_list,'False'))));
    
    %Calculate Z-based temperature and plot
    %Generate vector of temperatures for each point
    for idx=[1:length(time)]
        [Temp_seq(idx),Temp_er_seq(idx),~,~,~]=Z_temp_calculator_v5(spectral_radiance(idx,:)',wavelength',wl_exclude);
    end
    
    %Replace data with error larger than 1000K with NaN vales
    Temp_seq(Temp_er_seq>1000)=NaN;
    
    %Update tables
    handles.data_table.Data{file_idx,3}='True';
    handles.main_GUI_data.central_table.Data{file_idx,8}='True';
    
    %Update plots
    handles.rad_plot.XData= time;handles.rad_plot.YData = radiance;
    handles.temp_plot(2).XData= time;handles.temp_plot(2).YData = Temp_seq';
    handles.temp_plot(1).XData= time;handles.temp_plot(1).YData = zeros(1,length(time));
    
    %Update main data cell
    handles.main_GUI_data.main_data_table{file_idx,8}=Temp_seq';
    handles.main_GUI_data.main_data_table{file_idx,9}=Temp_er_seq';
    
end

guidata(hObject,handles);

% --------------------------------------------------------------------
function Z_all_Callback(hObject, eventdata, handles)
% hObject    handle to Z_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_short_p=handles.main_GUI_data.PMT_short_p;
file_list=handles.main_GUI_data.file_list;wavelength=handles.wavelength;

%Check data table
table_data=get(handles.data_table,'Data');

for idx=1:length(file_list);
    
    %Check for radiance data. If not present, skip this file
    if table_data{idx,2} == 'True'
        
        %Generate file name
        PMT_short_file=[file_list{idx} '_short.txt.'];
        
        %Update current shot
        handles.selected_file_txt.String=file_list{idx};
        
        %import necessary data
        time=handles.main_GUI_data.main_data_table{idx,5};
        radiance=handles.main_GUI_data.main_data_table{idx,6};
        spectral_radiance=handles.main_GUI_data.main_data_table{idx,7};
        
        %Turn off warning about negative data
        warning('off','MATLAB:Axes:NegativeDataInLogAxis');
        
        %Get list of wavelengths to exclude
        wl_list=get(handles.WL_table,'Data');
        wl_exclude=find(not(cellfun('isempty', strfind(wl_list,'False'))));
        
        %Calculate Z-based temperature and plot
        %Generate vector of temperatures for each point
        parfor Tidx=[1:length(time)]
            [Temp_seq(Tidx),Temp_er_seq(Tidx),~,~,~]=Z_temp_calculator_v5(spectral_radiance(Tidx,:)',wavelength',wl_exclude);
        end
        
        %Replace data with error larger than 1000K with NaN vales
        Temp_seq(Temp_er_seq>1000)=NaN;
        
        
        
        %Update tables
        handles.data_table.Data{idx,3}='True';
        handles.main_GUI_data.central_table.Data{idx,8}='True';
        
        handles.rad_plot.XData= time;handles.rad_plot.YData = radiance;
        handles.temp_plot(2).XData= time;handles.temp_plot(2).YData = Temp_seq';
        handles.temp_plot(1).XData= time;handles.temp_plot(1).YData = zeros(1,length(time));
        
        %Update main data cell
        handles.main_GUI_data.main_data_table{idx,8}=Temp_seq';
        handles.main_GUI_data.main_data_table{idx,9}=Temp_er_seq';
        
        %clear variables to prevent errors
        clear time radiance Temp_seq Temp_er_seq
        
        drawnow
    end
end

guidata(hObject,handles);


% --- Executes when user attempts to close emission_processing_GUI_bkg.
function emission_processing_GUI_bkg_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to emission_processing_GUI_bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.decision]=GUI_pass;

switch handles.decision
    case 1
        main_GUI=findobj('Tag','Emission_spectrometer_GUI');
        guidata(main_GUI,handles.main_GUI_data);
        delete(handles.emission_processing_GUI_bkg);
    case 2
        delete(handles.emission_processing_GUI_bkg);
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
delete(handles.emission_processing_GUI_bkg);
