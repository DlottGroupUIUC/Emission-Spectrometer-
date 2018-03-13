function varargout = Emission_spectrometer_GUI_v5_main(varargin)
% EMISSION_SPECTROMETER_GUI_V5_MAIN MATLAB code for Emission_spectrometer_GUI_v5_main.fig
%      EMISSION_SPECTROMETER_GUI_V5_MAIN, by itself, creates a new EMISSION_SPECTROMETER_GUI_V5_MAIN or raises the existing
%      singleton*.
%
%      H = EMISSION_SPECTROMETER_GUI_V5_MAIN returns the handle to a new EMISSION_SPECTROMETER_GUI_V5_MAIN or the handle to
%      the existing singleton*.
%
%      EMISSION_SPECTROMETER_GUI_V5_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMISSION_SPECTROMETER_GUI_V5_MAIN.M with the given input arguments.
%
%      EMISSION_SPECTROMETER_GUI_V5_MAIN('Property','Value',...) creates a new EMISSION_SPECTROMETER_GUI_V5_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Emission_spectrometer_GUI_v5_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Emission_spectrometer_GUI_v5_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Emission_spectrometer_GUI_v5_main

% Last Modified by GUIDE v2.5 13-Mar-2018 16:36:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Emission_spectrometer_GUI_v5_main_OpeningFcn, ...
    'gui_OutputFcn',  @Emission_spectrometer_GUI_v5_main_OutputFcn, ...
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


% --- Executes just before Emission_spectrometer_GUI_v5_main is made visible.
function Emission_spectrometer_GUI_v5_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Emission_spectrometer_GUI_v5_main (see VARARGIN)

% Choose default command line output for Emission_spectrometer_GUI_v5_main
handles.output = hObject;


%%%%Initialize main data cell.
%data as follows: 1. shot# 2. speed (km/s) 3. speed error (km/s) 4.impact time (ns)
%5. time (s) 6. radiance (W sr-1 m-2) 7. spectral radiance (W sr-1 m-3)
%8. Z temperature fits (K) 9. Z error (K) %10. graybody model 11. gray temperature (K)
%12. gray temperature error (K) 13. gray phi 14. gray phi error 15.
%calibration
handles.main_data_table=cell(1,15);


%Set Handles
guidata(hObject, handles);

% UIWAIT makes Emission_spectrometer_GUI_v5_main wait for user response (see UIRESUME)
% uiwait(handles.figure_main);


% --- Outputs from this function are returned to the command line.
function varargout = Emission_spectrometer_GUI_v5_main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function central_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to central_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%initialize table
set(hObject,'Data',cell(1,9));
set(hObject, 'ColumnName', {'shot #','speed (km/s)','error (km/s)', 'impact time (ns)' , 'PMT','PDV','radiance','Z (temp)','graybody'});
set(hObject, 'ColumnFormat',{'numeric','numeric','numeric','numeric','logical','logical','logical','logical','logical'});



% --- Executes when selected cell(s) is changed in central_table.
function central_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to central_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
try
    selected_cells=eventdata.Indices(1);
    
    file_idx=handles.file_idx;max_length=handles.max_length;
    file_list=handles.file_list;
    
    %Move to the file in the list which was selected
    file_idx=selected_cells;
    
    %Update the PDV/PMT file list
    PDV_file=[file_list{file_idx},'_ch3.txt'];
    PMT_file=[file_list{file_idx},'.tdms'];
    
    %Set display
    set(handles.file_txt,'String',file_list{file_idx});
    
    handles.file_idx=file_idx;handles.max_length=max_length;
    handles.PDV_file=PDV_file;
    handles.PMT_file=PMT_file;
    guidata(hObject,handles);
catch
end



% --------------------------------------------------------------------
function fwhm_all_Callback(hObject, eventdata, handles)
% hObject    handle to fwhm_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_short_p=handles.PMT_short_p;file_idx=handles.file_idx;
file_list=handles.file_list;

%Check data table
timing_data=get(handles.central_table,'Data');

for idx=file_idx:length(file_list)
    
    %Update current shot
    set(handles.plotted_file_txt,'String',file_list{idx});
    
    %Generate file name
    PMT_short_file=[file_list{idx} '_short.txt.'];
    
    %Import data
    data=importdata(fullfile(PMT_short_p,PMT_short_file));
    
    %Separate into time and radiance
    time=data(:,1);
    spectral_radiance=data(:,[2:end]);
    radiance=sum(spectral_radiance,2)*1e-9;
    
    %Turn off warning about negative data
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    
    %Plot fits
    %Check for an open rise time figure
    if isempty(findobj('Tag','FWHM fig'))
        handles.fwhm_fig=figure('Name','Full width half max of radiance',...
            'Tag','FWHM fig','NumberTitle','off','units','normalized','outerposition',[1 0 1 1]);
        
        %Change to the open rise time figure if it exists
    else
        figure(handles.fwhm_fig);
    end
    semilogx(time,radiance,'k');ylabel('radiance');xlabel('time (s)');
    title(file_list{file_idx},'FontSize',32);
    xlim([1e-9 1e-5]);
    
    if idx==file_idx
        [time_bounds ~]=ginput(2);
    end
    
    [~, l_bound] = min(abs(time-time_bounds(1))); %closest index for lower time bound
    [~, r_bound] = min(abs(time-time_bounds(2))); %closest index for higher time bound
    
    [max_value ~]=max(radiance(l_bound:r_bound));
    
    time_interp=linspace(time(l_bound),time(r_bound),5000);
    radiance_interp=spline(time,radiance,time_interp);
    
    low_idx=find(radiance_interp>=max_value*0.5,1,'first');
    high_idx=find(radiance_interp>=max_value*0.5,1,'last');
    
    fwhm=time_interp(high_idx)-time_interp(low_idx);
    
    figure(handles.fwhm_fig);
    logplot=semilogx(time,radiance,'dk',time_interp(low_idx),radiance_interp(low_idx),'or',time_interp(high_idx),...
        radiance_interp(high_idx),'or',time_interp,radiance_interp,'k');
    
    logplot(2).MarkerSize=30;
    logplot(3).MarkerSize=30;
    
    title(file_list{idx},'FontSize',32);
    xlim([1e-9 1e-5]);
    
    %Force update
    drawnow
    %Wait for user input
    pause
    
    running_list(idx,1)=str2num(file_list{idx});
    running_list(idx,2)=fwhm;
    
end

assignin('base','fwhm',running_list);

guidata(hObject,handles);


% --------------------------------------------------------------------
function rise_time_sgl_Callback(hObject, eventdata, handles)
% hObject    handle to rise_time_sgl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_short_p=handles.PMT_short_p;file_idx=handles.file_idx;
file_list=handles.file_list;

%Check data table
timing_data=get(handles.central_table,'Data');

%Update current shot
set(handles.plotted_file_txt,'String',file_list{file_idx});

%Generate file name
PMT_short_file=[file_list{file_idx} '_short.txt.'];

%Import data
data=importdata(fullfile(PMT_short_p,PMT_short_file));

%Separate into time and radiance
time=data(:,1);
spectral_radiance=data(:,[2:end]);
radiance=sum(spectral_radiance,2)*1e-9;

%Turn off warning about negative data
warning('off','MATLAB:Axes:NegativeDataInLogAxis');

%Plot fits
%Check for an open rise time figure
if isempty(findobj('Tag','rise time fig'))
    handles.rise_time_fig=figure('Name','Rise time of radiance (10% - 90%)',...
        'Tag','rise time fig','NumberTitle','off','units','normalized','outerposition',[1 0 1 1]);
    
    %Change to the open rise time figure if it exists
else
    figure(handles.rise_time_fig);
end
semilogx(time,radiance,'k');ylabel('radiance');xlabel('time (s)');
title(file_list{file_idx},'FontSize',32);
xlim([1e-9 1e-5]);

[time_bounds ~]=ginput(2);

[~, l_bound] = min(abs(time-time_bounds(1))); %closest index for lower time bound
[~, r_bound] = min(abs(time-time_bounds(2))); %closest index for higher time bound

[max_value ~]=max(radiance(l_bound:r_bound));

time_interp=linspace(time(l_bound),time(r_bound),5000);
radiance_interp=spline(time,radiance,time_interp);

low_idx=find(radiance_interp>=max_value*0.1,1,'first');
high_idx=find(radiance_interp>=max_value*0.9,1,'first');

rising_edge=time_interp(high_idx)-time_interp(low_idx);

figure(handles.rise_time_fig);
logplot=semilogx(time,radiance,'dk',time_interp(low_idx),radiance_interp(low_idx),'or',time_interp(high_idx),...
    radiance_interp(high_idx),'or',time_interp,radiance_interp,'k');
logplot(2).MarkerSize=20;
logplot(3).MarkerSize=20;
title([file_list{file_idx},' rise time=', num2str(rising_edge)],'FontSize',32);
xlim([1e-9 1e-5]);

guidata(hObject,handles);


% --------------------------------------------------------------------
function rise_time_all_Callback(hObject, eventdata, handles)
% hObject    handle to rise_time_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_short_p=handles.PMT_short_p;file_idx=handles.file_idx;
file_list=handles.file_list;

%Check data table
timing_data=get(handles.central_table,'Data');

for idx=file_idx:length(file_list)
    
    %Update current shot
    set(handles.plotted_file_txt,'String',file_list{idx});
    
    %Generate file name
    PMT_short_file=[file_list{idx} '_short.txt.'];
    
    %Import data
    data=importdata(fullfile(PMT_short_p,PMT_short_file));
    
    %Separate into time and radiance
    time=data(:,1);
    spectral_radiance=data(:,[2:end]);
    radiance=sum(spectral_radiance,2)*1e-9;
    
    %Turn off warning about negative data
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    
    %Plot fits
    %Check for an open rise time figure
    if isempty(findobj('Tag','rise time fig'))
        handles.rise_time_fig=figure('Name','Rise time of radiance (10% - 90%)',...
            'Tag','rise time fig','NumberTitle','off','units','normalized','outerposition',[1 0 1 1]);
        
        %Change to the open rise time figure if it exists
    else
        figure(handles.rise_time_fig);
    end
    semilogx(time,radiance,'k');ylabel('radiance');xlabel('time (s)');
    title(file_list{file_idx},'FontSize',32);
    xlim([1e-9 1e-5]);
    
    if idx==file_idx
        [time_bounds ~]=ginput(2);
    end
    
    [~, l_bound] = min(abs(time-time_bounds(1))); %closest index for lower time bound
    [~, r_bound] = min(abs(time-time_bounds(2))); %closest index for higher time bound
    
    [max_value ~]=max(radiance(l_bound:r_bound));
    
    time_interp=linspace(time(l_bound),time(r_bound),5000);
    radiance_interp=spline(time,radiance,time_interp);
    
    low_idx=find(radiance_interp>=max_value*0.1,1,'first');
    high_idx=find(radiance_interp>=max_value*0.9,1,'first');
    
    rising_edge=time_interp(high_idx)-time_interp(low_idx);
    
    figure(handles.rise_time_fig);
    logplot=semilogx(time,radiance,'dk',time_interp(low_idx),radiance_interp(low_idx),'or',time_interp(high_idx),...
        radiance_interp(high_idx),'or',time_interp,radiance_interp,'k');
    logplot(2).MarkerSize=20;
    logplot(3).MarkerSize=20;
    
    title(file_list{idx},'FontSize',32);
    xlim([1e-9 1e-5]);
    
    %Force update
    drawnow
    %Wait for user input
    pause
    
    running_list(idx,1)=str2num(file_list{idx});
    running_list(idx,2)=rising_edge;
    
end

assignin('base','rising_edges',running_list);

guidata(hObject,handles);



% --------------------------------------------------------------------
function fall_time_all_Callback(hObject, eventdata, handles)
% hObject    handle to fall_time_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PMT_short_p=handles.PMT_short_p;file_idx=handles.file_idx;
file_list=handles.file_list;

for idx=file_idx:length(file_list)
    
    %Update current shot
    set(handles.plotted_file_txt,'String',file_list{idx});
    
    %Generate file name
    PMT_short_file=[file_list{idx} '_short.txt.'];
    
    %Import data
    data=importdata(fullfile(PMT_short_p,PMT_short_file));
    
    %Separate into time and radiance
    time=data(:,1);
    spectral_radiance=data(:,[2:end]);
    radiance=sum(spectral_radiance,2)*1e-9;
    
    %Turn off warning about negative data
    warning('off','MATLAB:Axes:NegativeDataInLogAxis');
    
    %Plot fits
    %Check for an open rise time figure
    if isempty(findobj('Tag','fall time fig'))
        handles.fall_time_fig=figure('Name','Fall time of radiance (10% - 90%)',...
            'Tag','fall time fig','NumberTitle','off','units','normalized','outerposition',[1 0 1 1]);
        
        %Change to the open rise time figure if it exists
    else
        figure(handles.fall_time_fig);
    end
    semilogx(time,radiance,'k');ylabel('radiance');xlabel('time (s)');
    title(file_list{file_idx},'FontSize',32);
    xlim([1e-9 1e-5]);
    
    if idx==1
        [time_bounds ~]=ginput(2);
    end
    
    [~, l_bound] = min(abs(time-time_bounds(1))); %closest index for lower time bound
    [~, r_bound] = min(abs(time-time_bounds(2))); %closest index for higher time bound
    
    [max_value ~]=max(radiance(l_bound:r_bound));
    
    time_interp=linspace(time(l_bound),time(r_bound),5000);
    radiance_interp=spline(time,radiance,time_interp);
    
    low_idx=find(radiance_interp>=max_value*0.1,1,'last');
    high_idx=find(radiance_interp>=max_value*0.9,1,'last');
    
    falling_edge=time_interp(low_idx)-time_interp(high_idx);
    
    figure(handles.fall_time_fig);
    semilogx(time,radiance,'dk',time_interp(low_idx),radiance_interp(low_idx),'or',time_interp(high_idx),...
        radiance_interp(high_idx),'or',time_interp,radiance_interp,'k');
    title(file_list{idx},'FontSize',32);
    xlim([1e-9 1e-5]);
    
    %Force update
    drawnow
    %Wait for user input
    pause
    
    running_list(idx,1)=str2num(file_list{idx});
    running_list(idx,2)=falling_edge;
    
end

assignin('base','falling_edges',running_list);

guidata(hObject,handles);


% --------------------------------------------------------------------
function temp_analysis_GUI_btn_Callback(hObject, eventdata, handles)
% hObject    handle to temp_analysis_GUI_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

temperature_analysis_GUI;

% --------------------------------------------------------------------
function file_header_Callback(hObject, eventdata, handles)
% hObject    handle to file_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_header_Callback(hObject, eventdata, handles)
% hObject    handle to load_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_header_Callback(hObject, eventdata, handles)
% hObject    handle to save_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function raw_load_Callback(hObject, eventdata, handles)
% hObject    handle to raw_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
central_data=get(handles.central_table,'Data');



%If data table is populated, ask if you really want to load new
if min(cellfun('isempty',central_data))==1
    check_load='Yes';
else
    check_load=load_checker_v5;
end

switch check_load
    case 'Yes'
        clear main_data_table central_data
        %Load files to analyze
        [r,p]=uigetfile('*3.txt;*.tdms;*short.txt','multiselect','on');
        
        %Make cell if only 1 file is selected
        if ischar(r)==1
            disp(true)
            r={r};
        end
        
        %Initialize file index for scrolling through the different files
        file_idx=1;
        max_length=length(r);
        clear file_list
        
        %Check if parent folder is from PDV, PMT, or PMT shortened data
        %Check for PMT folder
        if strcmp(p(end-3:end),'PMT\')==1;
            
            %Assign PDV, PMT, PMT shortened, and the parent path
            PDV_p=[p(1:end-4) 'PDV\'];
            PMT_p=p;
            PMT_short_p=[p(1:end-4) 'PMT_short\'];
            parent_p=p(1:end-4);
            
            %Generate file list
            for idx=1:max_length
                file_names=cell2mat(r(idx));
                file_list{idx}=file_names(1:end-5);
            end
            
            %Define file names for the alternate files
            PMT_file=[file_list{file_idx},'.tdms'];
            PDV_file=[file_list{file_idx},'_ch3.txt'];
            
            %Check for PDV folder
        elseif strcmp(p(end-3:end),'PDV\')==1;
            
            %Assign PDV, PMT, PMT shortened, and the parent path
            PDV_p=p;
            PMT_p=[p(1:end-4) 'PMT\'];
            PMT_short_p=[p(1:end-4) 'PMT_short\'];
            parent_p=p(1:end-4);
            
            %Generate file list
            for idx=1:max_length;
                file_names=cell2mat(r(idx));
                file_list{idx}=file_names(1:end-8);
            end
            
            %Define file names for the alternate files
            PDV_file=[file_list{file_idx},'_ch3.txt'];
            PMT_file=[file_list{file_idx},'.tdms'];
            
            %Check for shortened PMT folder
        elseif strcmp(p(end-9:end),'PMT_short\')==1;
            
            %Assign PDV, PMT, PMT shortened, and the parent path
            PDV_p=[p(1:end-10) 'PDV\'];
            PMT_p=[p(1:end-10) 'PMT\'];
            PMT_short_p=p;
            parent_p=p(1:end-10);
            
            %Generate file list
            for idx=1:max_length
                file_names=cell2mat(r(idx));
                file_list{idx}=file_names(1:end-10);
            end
            
            %Define file names for the alternate files
            PMT_file=[file_list{file_idx},'.tdms'];
            PDV_file=[file_list{file_idx},'_ch3.txt'];
            
        end
        
        
        %Generate list of files, speeds, and delay times and check for files
        for idx=1:max_length
            central_data{idx,1}=file_list{idx};
            handles.main_data_table{idx,1}=file_list{idx};
            
            %Extend central table data out
            central_data{idx,7}='False';
            central_data{idx,8}='False';
            central_data{idx,9}='False';
            
            %Check for PMT file
            if exist(fullfile(PMT_p,[file_list{idx},'.tdms'])) == 2
                central_data{idx,5}='True';
            else
                central_data{idx,5}='False';
            end
            
            %Check for PDV file
            if exist(fullfile(PDV_p,[file_list{idx},'_ch3.txt'])) == 2
                central_data{idx,6}='True';
            else
                central_data{idx,6}='False';
            end
            
        end
        
        %Apply list to table
        set(handles.central_table,'Data',central_data);
        
        handles.file_idx=file_idx;handles.max_length=max_length;handles.r=r;
        handles.p=p;handles.timing_list=central_data;handles.file_list=file_list;
        handles.PDV_p=PDV_p;handles.PMT_p=PMT_p;handles.PDV_file=PDV_file;
        handles.PMT_file=PMT_file;handles.parent_p=parent_p;
        handles.PMT_short_p=PMT_short_p;handles.specT_emiss=cell(1,1);
        
    case 'No'
end

guidata(hObject,handles);

% --------------------------------------------------------------------
function processed_load_Callback(hObject, eventdata, handles)
% hObject    handle to processed_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Grab table data to check population
central_data=get(handles.central_table,'Data');

%If data table is populated, ask if you really want to load new
if min(cellfun('isempty',central_data))==1
    check_load='Yes';
else
    check_load=load_checker_v5;
end

if check_load=='Yes'
    %clear out past data
    clear clear central_data
    
    %Initialize file index for scrolling through the different files
    file_idx=1;
    
    %Load files to analyze
    [r,parent_p]=uigetfile('Select processed data file to load');
    data_initial=load(fullfile(parent_p,r));
    data_block=data_initial.data_block;
    
    %make sure r is long enough to read through data
    while length(r)<17
        r=strcat('spacer',r);
    end
    
    %Generate folder names
    PDV_p=[parent_p,'PDV\'];
    PMT_p=[parent_p,'PMT\'];
    PMT_short_p=[parent_p,'PMT_short\'];
    
    %%%identify file type and load accordingly%%%
    %For speed files
    if strcmp(r(end-13:end),'speed data.mat')==1
        
        %convert data to central table
        central_table=data_block;
        data_table=data_block;
        %Generate file list and other relevant table data
        for q=1:size(data_block,1)
            
            %generate file list
            file_list{q,1}=data_block{q,1};
            
            %Define file names for the alternate files
            PMT_file=[file_list{q},'.tdms'];
            PDV_file=[file_list{q},'_ch3.txt'];
            
            %Determine if PMT and PDV files exist
            if exist([PMT_p,PMT_file])==2
                central_table{q,5}='True';
            else
                central_table{q,5}='False';
            end
            
            if exist([PDV_p,PDV_file])==2
                central_table{q,6}='True';
            else
                central_table{q,6}='False';
            end
            
            %Extend central table data out
            central_table{q,9}='False';
            
        end
        
       
        %For radiance files
    elseif strcmp(r(end-16:end),'radiance data.mat')==1
        
        %convert data to central table
        central_table(:,1)=data_block(:,1);
        
        %Convert data to main table
        data_table(:,1)=data_block(:,1);
        data_table(:,[5,6,7])=data_block(:,[2,3,4]);
        
        %Generate file list and other relevant table data
        for q=1:size(data_block,1)
            
            %generate file list
            file_list{q,1}=data_block{q,1};
            
            %Define file names for the alternate files
            PMT_file=[file_list{q},'.tdms'];
            PDV_file=[file_list{q},'_ch3.txt'];
            
            %Determine if PMT and PDV files exist
            if exist([PMT_p,PMT_file])==2
                central_table{q,5}='True';
            else
                central_table{q,5}='False';
            end
            
            if exist([PDV_p,PDV_file])==2
                central_table{q,6}='True';
            else
                central_table{q,6}='False';
            end
            
            %Determine if radiance data exist
            if isempty(data_block{q,3}) == 0
                central_table{q,7}='True';
            else
                central_table{q,7}='False';
            end
            
            %Extend central table data out
            central_table{q,9}='False';
            
        end
    
        %For Z fits
    elseif strcmp(r(end-9:end),'Z data.mat')==1
        
        %convert data to central table
        central_table(:,1)=data_block(:,1);
        
        %Convert data to main table
        data_table(:,1)=data_block(:,1);
        data_table(:,[5,6,7,8,9])=data_block(:,[2,3,4,5,6]);
        
        %Generate file list and other relevant table data
        for q=1:size(data_block,1)
            
            %generate file list
            file_list{q,1}=data_block{q,1};
            
            %Define file names for the alternate files
            PMT_file=[file_list{q},'.tdms'];
            PDV_file=[file_list{q},'_ch3.txt'];
            
            %Determine if PMT and PDV files exist
            if exist([PMT_p,PMT_file])==2
                central_table{q,5}='True';
            else
                central_table{q,5}='False';
            end
            
            if exist([PDV_p,PDV_file])==2
                central_table{q,6}='True';
            else
                central_table{q,6}='False';
            end
            
            %Determine if radiance data exist
            if isempty(data_block{q,6}) == 0
                central_table{q,7}='True';
            else
                central_table{q,7}='False';
            end
            
            %Determine if Z data exist
            if isempty(data_block{q,3}) == 0
                central_table{q,8}='True';
            else
                central_table{q,8}='False';
            end
            
            %Determine if graybody data exist
            if isempty(data_block{q,3}) == 0
                central_table{q,9}='True';
            else
                central_table{q,9}='False';
            end
            
            %Extend central table data out
            central_table{q,9}='False';
            
        end
        
        %for graybody fits
    elseif strcmp(r(end-16:end),'graybody data.mat')==1
        
         %convert data to central table
        central_table(:,1)=data_block(:,1);
        
        %Convert data to main table
        data_table(:,1)=data_block(:,1);
        data_table(:,[5,6,7,8,10,11,12,13,14])=data_block(:,[2,3,4,5,6,7,8,9,10]);
        
        %Generate file list and other relevant table data
        for q=1:size(data_block,1)
            
            %generate file list
            file_list{q,1}=data_block{q,1};
            
            %Define file names for the alternate files
            PMT_file=[file_list{q},'.tdms'];
            PDV_file=[file_list{q},'_ch3.txt'];
            
            %Determine if PMT and PDV files exist
            if exist([PMT_p,PMT_file])==2
                central_table{q,5}='True';
            else
                central_table{q,5}='False';
            end
            
            if exist([PDV_p,PDV_file])==2
                central_table{q,6}='True';
            else
                central_table{q,6}='False';
            end
            
            %Determine if radiance data exist
            if isempty(data_block{q,6}) == 0
                central_table{q,7}='True';
            else
                central_table{q,7}='False';
            end
            
            %Determine if Z data exist
            if isempty(data_block{q,8}) == 0
                central_table{q,8}='True';
            else
                central_table{q,8}='False';
            end
            
            %Determine if graybody data exist
            if isempty(data_block{q,3}) == 0
                central_table{q,9}='True';
            else
                central_table{q,9}='False';
            end
            
        end
            
                    %for all data
    elseif strcmp(r(end-11:end),'all data.mat')==1
        
         %convert data to central table
        central_table(:,1)=data_block(:,1);
        central_table(:,2)=data_block(:,2);
        central_table(:,3)=data_block(:,3);
        central_table(:,4)=data_block(:,4);
        
        %Convert data to main table
        data_table=data_block;
                
        %Generate file list and other relevant table data
        for q=1:size(data_block,1)
            
            %generate file list
            file_list{q,1}=data_block{q,1};
            
            %Define file names for the alternate files
            PMT_file=[file_list{q},'.tdms'];
            PDV_file=[file_list{q},'_ch3.txt'];
            
            %Determine if PMT and PDV files exist
            if exist([PMT_p,PMT_file])==2
                central_table{q,5}='True';
            else
                central_table{q,5}='False';
            end
            
            if exist([PDV_p,PDV_file])==2
                central_table{q,6}='True';
            else
                central_table{q,6}='False';
            end
            
            %Determine if radiance data exist
            if isempty(data_block{q,6}) == 0
                central_table{q,7}='True';
            else
                central_table{q,7}='False';
            end
            
            %Determine if Z data exist
            if isempty(data_block{q,8}) == 0
                central_table{q,8}='True';
            else
                central_table{q,8}='False';
            end
            
            
            %Determine if graybody data exist
            if isempty(data_block{q,11}) == 0
                central_table{q,9}='True';
            else
                central_table{q,9}='False';
            end
            
           
        end
        end
        
               
    %Apply list to table
    handles.central_table.Data=central_table;
    handles.file_idx=file_idx;
    handles.p=parent_p;handles.file_list=file_list;
    handles.PDV_p=PDV_p;handles.PMT_p=PMT_p;handles.parent_p=parent_p;
    handles.PMT_short_p=PMT_short_p;handles.main_data_table=data_table;
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function Process_GUI_header_Callback(hObject, eventdata, handles)
% hObject    handle to Process_GUI_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function speed_timing_GUI_menu_Callback(hObject, eventdata, handles)
% hObject    handle to speed_timing_GUI_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
speed_timing_GUI;

% --------------------------------------------------------------------
function emission_GUI_menu_Callback(hObject, eventdata, handles)
% hObject    handle to emission_GUI_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
emission_processing_GUI;


% --------------------------------------------------------------------
function binning_GUI_menu_Callback(hObject, eventdata, handles)
% hObject    handle to binning_GUI_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
binning_GUI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 SAVING DATA AS both .MAT or .TXT                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function speed_save_Callback(hObject, eventdata, handles)
% hObject    handle to speed_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Get data to be saved
parent_p=handles.parent_p;file_list=handles.file_list;

data_block=handles.main_data_table(:,[1,2,3,4]);

%make sure data block is large enough to load
if size(data_block,2)<15
    data_block{1,15}=[];
end

save_flag=1;

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

text_name=[name,' speed data.txt'];

save_with_flags_v5(parent_p,file_list,data_block,text_name,save_flag);

mat_name=[name,' speed data.mat'];

save(fullfile(parent_p,mat_name),'data_block');

% --------------------------------------------------------------------
function radiance_save_Callback(hObject, eventdata, handles)
% hObject    handle to radiance_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parent_p=handles.parent_p;file_list=handles.file_list;

data_block=handles.main_data_table(:,[1,5,6,7]);

%make sure data block is large enough to load
if size(data_block,2)<15
    data_block{1,15}=[];
end

save_flag=2;

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

text_name=[name,' radiance data.txt'];

save_with_flags_v5(parent_p,file_list,data_block,text_name,save_flag);

mat_name=[name,' radiance data.mat'];

save(fullfile(parent_p,mat_name),'data_block');

% --------------------------------------------------------------------
function Z_save_Callback(hObject, eventdata, handles)
% hObject    handle to Z_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parent_p=handles.parent_p;file_list=handles.file_list;

data_block=handles.main_data_table(:,[1,5,6,7,8,9]);

%make sure data block is large enough to load
if size(data_block,2)<15
    data_block{1,15}=[];
end

save_flag=3;

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

text_name=[name,' Z data.txt'];

save_with_flags_v5(parent_p,file_list,data_block,text_name,save_flag);

mat_name=[name,' Z data.mat'];

save(fullfile(parent_p,mat_name),'data_block');

% --------------------------------------------------------------------
function gray_save_Callback(hObject, eventdata, handles)
% hObject    handle to gray_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parent_p=handles.parent_p;file_list=handles.file_list;

data_block=handles.main_data_table(:,[1,5,6,7,8,10,11,12,13,14]);

%make sure data block is large enough to load
if size(data_block,2)<15
    data_block{1,15}=[];
end

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

%Save just temperatures
save_flag=4;
text_T_name=[name,' gray temperature data.txt'];

save_with_flags_v5(parent_p,file_list,data_block,text_T_name,save_flag);

%Save just phi
save_flag=5;
text_T_name=[name,' gray phi data.txt'];

save_with_flags_v5(parent_p,file_list,data_block,text_T_name,save_flag);

mat_name=[name,' graybody data.mat'];

save(fullfile(parent_p,mat_name),'data_block');

% --------------------------------------------------------------------
function all_save_Callback(hObject, eventdata, handles)
% hObject    handle to all_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parent_p=handles.parent_p;file_list=handles.file_list;

data_block=handles.main_data_table;

%make sure data block is large enough to load
if size(data_block,2)<15
    data_block{1,15}=[];
end

%Prompt for file name
prompt={'Enter file save name'};
dlg_title='Saved file name';
name=[cell2mat(inputdlg(prompt,dlg_title))];

%Save only matlab file for all save
mat_name=[name,' all data.mat'];

save(fullfile(parent_p,mat_name),'data_block');


% --------------------------------------------------------------------
function analyze_header_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function analyze_gray_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
graybody_spectrum_analysis_GUI;
