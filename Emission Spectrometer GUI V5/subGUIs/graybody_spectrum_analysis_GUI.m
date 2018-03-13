function varargout = graybody_spectrum_analysis_GUI(varargin)
% GRAYBODY_SPECTRUM_ANALYSIS_GUI MATLAB code for graybody_spectrum_analysis_GUI.fig
%      GRAYBODY_SPECTRUM_ANALYSIS_GUI, by itself, creates a new GRAYBODY_SPECTRUM_ANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = GRAYBODY_SPECTRUM_ANALYSIS_GUI returns the handle to a new GRAYBODY_SPECTRUM_ANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      GRAYBODY_SPECTRUM_ANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRAYBODY_SPECTRUM_ANALYSIS_GUI.M with the given input arguments.
%
%      GRAYBODY_SPECTRUM_ANALYSIS_GUI('Property','Value',...) creates a new GRAYBODY_SPECTRUM_ANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before graybody_spectrum_analysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to graybody_spectrum_analysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help graybody_spectrum_analysis_GUI

% Last Modified by GUIDE v2.5 13-Mar-2018 17:16:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @graybody_spectrum_analysis_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @graybody_spectrum_analysis_GUI_OutputFcn, ...
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


% --- Executes just before graybody_spectrum_analysis_GUI is made visible.
function graybody_spectrum_analysis_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to graybody_spectrum_analysis_GUI (see VARARGIN)

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
handles.data_table.Data=main_GUI_data.central_table.Data(:,[1,7,9]);

%Assign center wavelength for fiber and calculate inverse wavelengths
wavelength=[444.88,450.88,457.38,464.38,471.88,480,488.88,498.5,509,520.5,...
    530.88,540,549.75,560.13,571.38,583.63,596.88,611.25,623.75,634,644.88,...
    656.5,668.88,682,696,711,727.13,744.5,763.25,783.5,805.38,829.25]';

%Set current plot
handles.selected_file_txt.String=main_GUI_data.file_list{file_idx};

%%%Initialize Plots%%%
%Generate data
%import necessary data
if sum(cellfun('isempty',main_GUI_data.main_data_table(1,[5,6,7,8,10,11,13])))==0;
    time=main_GUI_data.main_data_table{file_idx,5};
    radiance=main_GUI_data.main_data_table{file_idx,6};
    spectral_radiance=main_GUI_data.main_data_table{file_idx,7}(25,:);
    T=main_GUI_data.main_data_table{file_idx,11};
    Ter=main_GUI_data.main_data_table{file_idx,12};
    P=main_GUI_data.main_data_table{file_idx,13};
    Per=main_GUI_data.main_data_table{file_idx,14};
    model_SR=main_GUI_data.main_data_table{file_idx,10}(25,:);
    
else
    time=zeros(1,1e2);
    radiance=zeros(1,1e2);
    spectral_radiance=zeros(1,32);
    T=zeros(1e2);
    P=zeros(1e2);
    model_SR=zeros(1,32);
end

%initialize plots with real data
figure(handles.graybody_spectrum_analysis_GUI_bkg);
handles.WL_sub=subplot('Position',[0.20,0.63 0.65 0.35]);handles.WL_plot=plot(wavelength,spectral_radiance,'or',wavelength,model_SR,'-k');
handles.TP_sub=subplot('Position',[0.20 0.25 0.65 0.30]);
yyaxis left; handles.TP_plot(1)=semilogx(time,T,'ok');ylabel('temperature (K)');ylim([0 7000]);
yyaxis right; handles.TP_plot(2)=semilogx(time,P,'xr');ylabel('\phi');
handles.rad_sub=subplot('Position',[0.20 0.15 0.65 0.10]);handles.rad_plot=semilogx(time,radiance,'b',time(25),radiance(25),'dr');

try
legend(handles.WL_sub,['t=',num2str(round(time(25)*1e9,1)),' ns'],...
    ['T=',num2str(round(T(25)/100)*100),'±',num2str(round(Ter(25)/10)*10),'K',...
    '\epsilon=',num2str(round(P(25)*1000)/1000),'±',num2str(round(Per(25)*1000)/1000)],...
    'Location','Northwest');
catch
end

%adjust plot windows
set(handles.TP_sub,'XTickLabel',[]);
ylabel(handles.rad_sub,'radiance (W/sr \cdot m^2)');
ylabel(handles.WL_sub,'spectral radiance (W/sr \cdot m^3)');
xlabel(handles.WL_sub,'wavelength (nm)');
xlabel(handles.rad_sub,'time (s)');
linkaxes([handles.rad_sub handles.TP_sub],'x');
xlim(handles.rad_sub,[1e-9 1e-5]);



%initialize autofit flag
handles.auto_flag.Value=1;

%Set callback for auto fitting
handles.TP_sub.ButtonDownFcn=@(hObject,eventdata)graybody_spectrum_analysis_GUI('TP_sub_ButtonDownFcn',hObject,eventdata,guidata(hObject));
handles.rad_sub.ButtonDownFcn=@(hObject,eventdata)graybody_spectrum_analysis_GUI('rad_sub_ButtonDownFcn',hObject,eventdata,guidata(hObject));

%Set hittest to off so clicking on lines is ok
handles.TP_plot(1).HitTest='off';
handles.TP_plot(2).HitTest='off';
handles.rad_plot(1).HitTest='off';
handles.rad_plot(2).HitTest='off';

%initialize close function
handles.decision=3;

% Update main GUI
handles.output = hObject;
handles.main_GUI_data=main_GUI_data;
handles.file_idx=file_idx;
handles.PMT_short_file=PMT_short_file;
handles.wavelength=wavelength;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes graybody_spectrum_analysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.graybody_spectrum_analysis_GUI_bkg);


% --- Outputs from this function are returned to the command line.
function varargout = graybody_spectrum_analysis_GUI_OutputFcn(hObject, eventdata, handles)
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
%initialize table
set(hObject,'Data',cell(1,3));
set(hObject, 'ColumnName', {'shot #','radiance','graybody' });
set(hObject, 'ColumnFormat',{'numeric','logical','logical',});


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
    
    %Check if autofit is on
    if handles.auto_flag.Value==1
        
        %retreive data for replotting
        if sum(cellfun('isempty',handles.main_GUI_data.main_data_table(1,[5,6,7,10,11,13])))==0;
            time=handles.main_GUI_data.main_data_table{file_idx,5};
            radiance=handles.main_GUI_data.main_data_table{file_idx,6};
            spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7}(25,:);
            T=handles.main_GUI_data.main_data_table{file_idx,11};
            Ter=handles.main_GUI_data.main_data_table{file_idx,12};
            P=handles.main_GUI_data.main_data_table{file_idx,13};
            Per=handles.main_GUI_data.main_data_table{file_idx,14};
            model_SR=handles.main_GUI_data.main_data_table{file_idx,10}(25,:);
            
        else
            time=zeros(1,1e2);
            radiance=zeros(1,1e2);
            spectral_radiance=zeros(1,32);
            T=zeros(1e2);
            P=zeros(1e2);
        end
        
        %update plots
        handles.rad_plot(1).XData= time;handles.rad_plot(1).YData = radiance;
        handles.rad_plot(2).XData=time(25);handles.rad_plot(2).YData=radiance(25);
        handles.TP_plot(1).XData= time;handles.TP_plot(1).YData = T;
        handles.TP_plot(2).XData= time;handles.TP_plot(2).YData = P;
        handles.WL_plot(1).YData=spectral_radiance;handles.WL_plot(2).YData=model_SR;
        
        try
            legend(handles.WL_sub,['t=',num2str(round(time(25)*1e9,1)),' ns'],...
                ['T=',num2str(round(T(25)/100)*100),'±',num2str(round(Ter(25)/10)*10),'K',...
                '\epsilon=',num2str(round(P(25)*1000)/1000),'±',num2str(round(Per(25)*1000)/1000)],...
                'Location','Northwest');
        catch
        end
        
        
    end
    
    handles.file_idx=file_idx;
    handles.PMT_short_file=PMT_short_file;
    guidata(hObject,handles);
catch
end


% --- Executes when user attempts to close graybody_spectrum_analysis_GUI_bkg.
function graybody_spectrum_analysis_GUI_bkg_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to graybody_spectrum_analysis_GUI_bkg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.decision]=GUI_pass;

switch handles.decision
    case 1
        main_GUI=findobj('Tag','Emission_spectrometer_GUI');
        guidata(main_GUI,handles.main_GUI_data);
        delete(handles.graybody_spectrum_analysis_GUI_bkg);
    case 2
        delete(handles.graybody_spectrum_analysis_GUI_bkg);
    case 3
end


% --------------------------------------------------------------------
function save_close_Callback(hObject, eventdata, handles)
% hObject    handle to save_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Grab fits
stored_fits=handles.fits_table.Data;
file_list=handles.main_GUI_data.file_list;

%Convert fit numbers to numerics to aid sorting
stored_fits(:,1)=cellfun(@str2num,stored_fits(:,1),'UniformOutput',0);
fit_idx=cellfun('isempty',stored_fits(:,1));

%Convert non-numeric names to NaN if necessary. They will be sorted to the
%bottom
if sum(fit_idx)>=1
stored_fits{fit_idx,1}=NaN;
end

%Sort fits by shot number and time
sorted_fits=sortrows(stored_fits,[1,2]);
sorted_fits(:,1)=cellfun(@num2str,sorted_fits(:,1),'UniformOutput',0);

%grab wavelength and initialize cell for saving
fit_save_cell{1,1}=handles.wavelength';

%Gather spectra for saving
for idx=1:size(sorted_fits,1);
    
    shot=sorted_fits{idx,1};
    time_pt=sorted_fits{idx,2}*1e-9;
    
    if shot~= NaN
        
        %Move to the file in the list which was selected
        idxC=strfind(file_list,shot);
        file_idx = find(not(cellfun('isempty', idxC)));
        
        %retreive data for storing
        time=handles.main_GUI_data.main_data_table{file_idx,5};
        
        %locate time index
        time_idx=find(time>=time_pt,1,'first');
        
        spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7}(time_idx,:);
        T=handles.main_GUI_data.main_data_table{file_idx,11}(time_idx);
        Ter=handles.main_GUI_data.main_data_table{file_idx,12}(time_idx);
        P=handles.main_GUI_data.main_data_table{file_idx,13}(time_idx);
        Per=handles.main_GUI_data.main_data_table{file_idx,14}(time_idx);
        model_SR=handles.main_GUI_data.main_data_table{file_idx,10}(time_idx,:);
        
        %Define cell for fitting
        fit_save_cell{idx+1,1}=spectral_radiance;
        fit_save_cell{idx+1,2}=model_SR;
        fit_save_cell{idx+1,3}=[T,Ter,P,Per];
        fit_save_cell{idx+1,4}=['shot ',shot,' @ ',num2str(round(time_pt*1e9,1)),'ns'];
        
    end
    
end

parent_p=handles.main_GUI_data.parent_p;

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

text_name=[name,' fit set.txt'];

save_T_fits_v5(parent_p,fit_save_cell,text_name);

guidata(hObject,handles);

main_GUI=findobj('Tag','Emission_spectrometer_GUI');
guidata(main_GUI,handles.main_GUI_data);
delete(handles.graybody_spectrum_analysis_GUI_bkg);




% --- Executes on button press in auto_flag.
function auto_flag_Callback(hObject, eventdata, handles)
% hObject    handle to auto_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_flag


% --------------------------------------------------------------------
function fit_ns_peak_header_Callback(hObject, eventdata, handles)
% hObject    handle to fit_ns_peak_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function store_fit_Callback(hObject, eventdata, handles)
% hObject    handle to store_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        %Find current size of the table
        
        
        if isempty(handles.fits_table.Data{1,1})==1
            table_size=0;
        else
            table_size=size(handles.fits_table.Data,1);
        end
        
        %Add shot record to current size
        handles.fits_table.Data{table_size+1,1}=handles.main_GUI_data.file_list{handles.file_idx};
        handles.fits_table.Data{table_size+1,2}=handles.rad_plot(2).XData*1e9; 


% --------------------------------------------------------------------
function fit_ns_sgl_Callback(hObject, eventdata, handles)
% hObject    handle to fit_ns_sgl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;
file_list=handles.main_GUI_data.file_list;

%retreive data for replotting
if sum(cellfun('isempty',handles.main_GUI_data.main_data_table(file_idx,[5,6,7,10,11,13])))==0;
    time=handles.main_GUI_data.main_data_table{file_idx,5};
    radiance=handles.main_GUI_data.main_data_table{file_idx,6};

    
else
    time=zeros(1,1e2);
    radiance=zeros(1,1e2);
    spectral_radiance=zeros(1,32);
    T=zeros(1e2);
    P=zeros(1e2);
end


%Copy data to manipulate
delay_line=radiance;
delay_time=time;

%Clip to 100 nanoseconds
delay_line(delay_time>1e-7)=[];
delay_time(delay_time>1e-7)=[];


peak_min=delay_line(1);

%Turn off warning about invalid min peak height
warning('off','signal:findpeaks:largeMinPeakHeight');

[~,locs] = findpeaks(delay_line,'MinPeakHeight',2*peak_min,...
    'MinPeakDistance',5e-9);

%Check for a found peak
if isempty(locs)==0;
    
    [~, loc_max]=max(delay_line(locs));
    
    xidx=locs(loc_max);
    
    spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7}(xidx,:);
    T=handles.main_GUI_data.main_data_table{file_idx,11};
    Ter=handles.main_GUI_data.main_data_table{file_idx,12};
    P=handles.main_GUI_data.main_data_table{file_idx,13};
    Per=handles.main_GUI_data.main_data_table{file_idx,14};
    model_SR=handles.main_GUI_data.main_data_table{file_idx,10}(xidx,:);
    
    %update plots
    handles.rad_plot(1).XData= time;handles.rad_plot(1).YData = radiance;
    handles.rad_plot(2).XData= time(xidx);handles.rad_plot(2).YData = radiance(xidx);
    handles.TP_plot(1).XData= time;handles.TP_plot(1).YData = T;
    handles.TP_plot(2).XData= time;handles.TP_plot(2).YData = P;
    handles.WL_plot(1).YData=spectral_radiance;handles.WL_plot(2).YData=model_SR;
    

    legend(handles.WL_sub,['t=',num2str(round(time(xidx)*1e9,1)),' ns'],...
        ['T=',num2str(round(T(xidx)/100)*100),'±',num2str(round(Ter(xidx)/10)*10),'K',...
        '\epsilon=',num2str(round(P(xidx)*1000)/1000),'±',num2str(round(Per(xidx)*1000)/1000)],...
        'Location','Northwest');


    
    %Find current size of the table
    if isempty(handles.fits_table.Data{1,1})==1
        table_size=0;
    else
        table_size=size(handles.fits_table.Data,1);
    end
    
    %Add shot record to current size
    handles.fits_table.Data{table_size+1,1}=handles.main_GUI_data.file_list{file_idx};
    handles.fits_table.Data{table_size+1,2}=time(xidx)*1e9;
      
    
end
    
    
guidata(hObject,handles);

% --------------------------------------------------------------------
function fit_ns_all_Callback(hObject, eventdata, handles)
% hObject    handle to fit_ns_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_list=handles.main_GUI_data.file_list;

for idx=1:length(file_list)
    
    %retreive data for replotting
    if sum(cellfun('isempty',handles.main_GUI_data.main_data_table(idx,[5,6,7,10,11,13])))==0;
        time=handles.main_GUI_data.main_data_table{idx,5};
        radiance=handles.main_GUI_data.main_data_table{idx,6};
        
        
    else
        time=zeros(1,1e2);
        radiance=zeros(1,1e2);
        spectral_radiance=zeros(1,32);
        T=zeros(1e2);
        P=zeros(1e2);
    end
    
    
    %Copy data to manipulate
    delay_line=radiance;
    delay_time=time;
    
    %Clip to 100 nanoseconds
    delay_line(delay_time>1e-7)=[];
    delay_time(delay_time>1e-7)=[];
    
    
    peak_min=delay_line(1);
    
    %Turn off warning about invalid min peak height
    warning('off','signal:findpeaks:largeMinPeakHeight');
    
    [~,locs] = findpeaks(delay_line,'MinPeakHeight',2*peak_min,...
        'MinPeakDistance',5e-9);
    
    %Check for a found peak
    if isempty(locs)==0;
        
        [~, loc_max]=max(delay_line(locs));
        
        xidx=locs(loc_max);
        
        spectral_radiance=handles.main_GUI_data.main_data_table{idx,7}(xidx,:);
        T=handles.main_GUI_data.main_data_table{idx,11};
        Ter=handles.main_GUI_data.main_data_table{idx,12};
        P=handles.main_GUI_data.main_data_table{idx,13};
        Per=handles.main_GUI_data.main_data_table{idx,14};
        model_SR=handles.main_GUI_data.main_data_table{idx,10}(xidx,:);
        
        %update plots
        handles.rad_plot(1).XData= time;handles.rad_plot(1).YData = radiance;
        handles.rad_plot(2).XData= time(xidx);handles.rad_plot(2).YData = radiance(xidx);
        handles.TP_plot(1).XData= time;handles.TP_plot(1).YData = T;
        handles.TP_plot(2).XData= time;handles.TP_plot(2).YData = P;
        handles.WL_plot(1).YData=spectral_radiance;handles.WL_plot(2).YData=model_SR;
        
        
        legend(handles.WL_sub,['t=',num2str(round(time(xidx)*1e9,1)),' ns'],...
            ['T=',num2str(round(T(xidx)/100)*100),'±',num2str(round(Ter(xidx)/10)*10),'K',...
            '\epsilon=',num2str(round(P(xidx)*1000)/1000),'±',num2str(round(Per(xidx)*1000)/1000)],...
            'Location','Northwest');
        
        
        
        %Find current size of the table
        if isempty(handles.fits_table.Data{1,1})==1
            table_size=0;
        else
            table_size=size(handles.fits_table.Data,1);
        end
        
        %Add shot record to current size
        handles.fits_table.Data{table_size+1,1}=handles.main_GUI_data.file_list{idx};
        handles.fits_table.Data{table_size+1,2}=time(xidx)*1e9;
        
    end
    
    drawnow
    pause(0.5);
    
end


% --- Executes on mouse press over axes background.
function TP_sub_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to temp_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if autofit is on
if handles.auto_flag.Value==1
    handles = guidata( ancestor(hObject, 'figure') );
    
    file_idx=handles.file_idx;
    
    
    %Identify cursor location
    xpos=eventdata.IntersectionPoint(1);
    ypos=eventdata.IntersectionPoint(2);
     
    %Retrieve data for replotting
    time=handles.main_GUI_data.main_data_table{file_idx,5};
    
    %Update Graybody fit
    %locate time index
    time_idx=find(time<xpos,1,'last');
    
    radiance=handles.main_GUI_data.main_data_table{file_idx,6};
    spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7}(time_idx,:);
    T=handles.main_GUI_data.main_data_table{file_idx,11};
    Ter=handles.main_GUI_data.main_data_table{file_idx,12};
    P=handles.main_GUI_data.main_data_table{file_idx,13};
    Per=handles.main_GUI_data.main_data_table{file_idx,14};
    model_SR=handles.main_GUI_data.main_data_table{file_idx,10}(time_idx,:);
    
    %update plots
    handles.rad_plot(1).XData= time;handles.rad_plot(1).YData = radiance;
    handles.rad_plot(2).XData= time(time_idx);handles.rad_plot(2).YData = radiance(time_idx);
    handles.TP_plot(1).XData= time;handles.TP_plot(1).YData = T;
    handles.TP_plot(2).XData= time;handles.TP_plot(2).YData = P;
    handles.WL_plot(1).YData=spectral_radiance;handles.WL_plot(2).YData=model_SR;
    
    legend(handles.WL_sub,['t=',num2str(round(time(time_idx)*1e9,1)),' ns'],...
        ['T=',num2str(round(T(time_idx)/100)*100),'±',num2str(round(Ter(time_idx)/10)*10),'K',...
        '\epsilon=',num2str(round(P(time_idx)*1000)/1000),'±',num2str(round(Per(time_idx)*1000)/1000)],...
        'Location','Northwest');
    

    
end

guidata(hObject,handles);

% --- Executes on mouse press over axes background.
function rad_sub_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to temp_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TP_sub_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fits_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fits_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% %initialize table
set(hObject,'Data',cell(1,2));
set(hObject, 'ColumnName', {'shot','time (ns)' });
set(hObject, 'ColumnFormat',{'numeric','numeric',});


% --- Executes when selected cell(s) is changed in fits_table.
function fits_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to fits_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
file_list=handles.main_GUI_data.file_list;

try
    shot=handles.fits_table.Data{eventdata.Indices(1),1};
    time_pt=handles.fits_table.Data{eventdata.Indices(1),2}*1e-9;
    
    %Move to the file in the list which was selected
    idxC=strfind(file_list,shot);
    file_idx = find(not(cellfun('isempty', idxC)));
        
    %Set display
    set(handles.selected_file_txt,'String',file_list{file_idx});
    
    %Check if autofit is on
    if handles.auto_flag.Value==1
        
        %retreive data for replotting
        time=handles.main_GUI_data.main_data_table{file_idx,5};
        
        %locate time index
        time_idx=find(time>=time_pt,1,'first');
        
        radiance=handles.main_GUI_data.main_data_table{file_idx,6};
        spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7}(time_idx,:);
        T=handles.main_GUI_data.main_data_table{file_idx,11};
        Ter=handles.main_GUI_data.main_data_table{file_idx,12};
        P=handles.main_GUI_data.main_data_table{file_idx,13};
        Per=handles.main_GUI_data.main_data_table{file_idx,14};
        model_SR=handles.main_GUI_data.main_data_table{file_idx,10}(time_idx,:);
        
        
        
        %update plots
        handles.rad_plot(1).XData= time;handles.rad_plot(1).YData = radiance;
        handles.rad_plot(2).XData=time(time_idx);handles.rad_plot(2).YData=radiance(time_idx);
        handles.TP_plot(1).XData= time;handles.TP_plot(1).YData = T;
        handles.TP_plot(2).XData= time;handles.TP_plot(2).YData = P;
        handles.WL_plot(1).YData=spectral_radiance;handles.WL_plot(2).YData=model_SR;
        
        try
            legend(handles.WL_sub,['t=',num2str(round(time(time_idx)*1e9,1)),' ns'],...
                ['T=',num2str(round(T(time_idx)/100)*100),'±',num2str(round(Ter(time_idx)/10)*10),'K',...
                '\epsilon=',num2str(round(P(time_idx)*1000)/1000),'±',num2str(round(Per(time_idx)*1000)/1000)],...
                'Location','Northwest');
        catch
        end
        
        
    end
    
    %Keep selected cells in case of deletion
%     if size
    handles.stored_selection=eventdata.Indices(:,1);
    handles.file_idx=file_idx;
    guidata(hObject,handles);
catch
end


% --------------------------------------------------------------------
function delete_shots_header_Callback(hObject, eventdata, handles)
% hObject    handle to delete_shots_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dlt_shots_select_Callback(hObject, eventdata, handles)
% hObject    handle to dlt_shots_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)  
    handles.fits_table.Data(handles.stored_selection,:)=[];
    
    
% --------------------------------------------------------------------
function dlt_shots_all_Callback(hObject, eventdata, handles)
% hObject    handle to dlt_shots_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fits_table.Data=cell(1,2);


% --------------------------------------------------------------------
function recalc_Callback(hObject, eventdata, handles)
% hObject    handle to recalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_idx=handles.file_idx;
file_list=handles.main_GUI_data.file_list;wavelength=handles.wavelength;


%Check data table
table_data=get(handles.data_table,'Data');

%check for graybody data. If not present, skip
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
    
    %set arbitrary first time to plot
    time_idx=25;
    
        
        %update plots
        handles.rad_plot(1).XData= time;handles.rad_plot(1).YData = radiance;
        handles.rad_plot(2).XData=time(time_idx);handles.rad_plot(2).YData=radiance(time_idx);
        handles.TP_plot(1).XData= time;handles.TP_plot(1).YData = spec_T;
        handles.TP_plot(2).XData= time;handles.TP_plot(2).YData = spec_emiss;
        handles.WL_plot(1).YData=spectral_radiance(time_idx,:);handles.WL_plot(2).YData=model(time_idx,:);
        
        try
            legend(handles.WL_sub,['t=',num2str(round(time(time_idx)*1e9,1)),' ns'],...
                ['T=',num2str(round(spec_T(time_idx)/100)*100),'±',num2str(round(spec_T_unc(time_idx)/10)*10),'K',...
                '\epsilon=',num2str(round(spec_emiss(time_idx)*1000)/1000),'±',num2str(round(spec_emiss_unc(time_idx)*1000)/1000)],...
                'Location','Northwest');
        catch
        end
    
    %Update main data cell
    handles.main_GUI_data.main_data_table{file_idx,10}=model;
    handles.main_GUI_data.main_data_table{file_idx,11}=spec_T';
    handles.main_GUI_data.main_data_table{file_idx,12}=spec_T_unc';
    handles.main_GUI_data.main_data_table{file_idx,13}=spec_emiss';
    handles.main_GUI_data.main_data_table{file_idx,14}=spec_emiss_unc';
    
    
end


guidata(hObject,handles);


% --------------------------------------------------------------------
function save_no_close_Callback(hObject, eventdata, handles)
% hObject    handle to save_no_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Grab fits
stored_fits=handles.fits_table.Data;
file_list=handles.main_GUI_data.file_list;

%Convert fit numbers to numerics to aid sorting
stored_fits(:,1)=cellfun(@str2num,stored_fits(:,1),'UniformOutput',0);
fit_idx=cellfun('isempty',stored_fits(:,1));

%Convert non-numeric names to NaN if necessary. They will be sorted to the
%bottom
if sum(fit_idx)>=1
stored_fits{fit_idx,1}=NaN;
end

%Sort fits by shot number and time
sorted_fits=sortrows(stored_fits,[1,2]);
sorted_fits(:,1)=cellfun(@num2str,sorted_fits(:,1),'UniformOutput',0);

%grab wavelength and initialize cell for saving
fit_save_cell{1,1}=handles.wavelength';

%Gather spectra for saving
for idx=1:size(sorted_fits,1);
    
    shot=sorted_fits{idx,1};
    time_pt=sorted_fits{idx,2}*1e-9;
    
    if shot~= NaN
        
        %Move to the file in the list which was selected
        idxC=strfind(file_list,shot);
        file_idx = find(not(cellfun('isempty', idxC)));
        
        %retreive data for storing
        time=handles.main_GUI_data.main_data_table{file_idx,5};
        
        %locate time index
        time_idx=find(time>=time_pt,1,'first');
        
        spectral_radiance=handles.main_GUI_data.main_data_table{file_idx,7}(time_idx,:);
        T=handles.main_GUI_data.main_data_table{file_idx,11}(time_idx);
        Ter=handles.main_GUI_data.main_data_table{file_idx,12}(time_idx);
        P=handles.main_GUI_data.main_data_table{file_idx,13}(time_idx);
        Per=handles.main_GUI_data.main_data_table{file_idx,14}(time_idx);
        model_SR=handles.main_GUI_data.main_data_table{file_idx,10}(time_idx,:);
        
        %Define cell for fitting
        fit_save_cell{idx+1,1}=spectral_radiance;
        fit_save_cell{idx+1,2}=model_SR;
        fit_save_cell{idx+1,3}=[T,Ter,P,Per];
        fit_save_cell{idx+1,4}=['shot ',shot,' @ ',num2str(round(time_pt*1e9,1)),'ns'];
        
    end
    
end

parent_p=handles.main_GUI_data.parent_p;

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

text_name=[name,' fit set.txt'];

save_T_fits_v5(parent_p,fit_save_cell,text_name);

guidata(hObject,handles);


% --------------------------------------------------------------------
function pass_close_Callback(hObject, eventdata, handles)
% hObject    handle to pass_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_GUI=findobj('Tag','Emission_spectrometer_GUI');
guidata(main_GUI,handles.main_GUI_data);
delete(handles.graybody_spectrum_analysis_GUI_bkg);
