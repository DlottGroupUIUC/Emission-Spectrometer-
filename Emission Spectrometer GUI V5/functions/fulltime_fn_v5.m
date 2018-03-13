function [out_t,out_speed,out_time,out_intensity] = fulltime_fn_v4(file_name,path,channel)
%fulltime_fn A functional form of the fulltime spline script for use with
%GUIs

%set resolution
res=160;
%Set minimum amplitude modulation for peak finding
s=.014;
 
%Select which channel to load
switch channel
    case 1
        %load channel 1 to analyze
        file=importdata(fullfile(path,[file_name,'_ch1.txt']));
        time=file.data(:,1)+abs(file.data(1,1)); %shifts beginning time to 0
        amp=file.data(:,2);
        
    case 2
        %load channel 2 to analyze
        file=importdata(fullfile(path,[file_name,'_ch2.txt']));
        time=file.data(:,1)+abs(file.data(1,1)); %shifts beginning time to 0
        amp=file.data(:,2);
        
    case 4
        %load channel 4 to analyze
        file=importdata(fullfile(path,[file_name,'_ch4.txt']));
        time=file.data(:,1)+abs(file.data(1,1)); %shifts beginning time to 0
        amp=file.data(:,2);
end

% Determine Sampling Frequency and Calculate Resolution
fs=1./(time(2));
r=floor((fs./(res.*10^6))./2); %sets r as the number of points needed
if mod(r,2)
else r=r-1;
end

%Load Trigger Data on Channel 3
file2=importdata(fullfile(path,[file_name,'_ch3.txt']));
trigt=time;  %separate data file into time
trigs=file2.data(:,2);  %and signal variables
sden=wden(trigs,'heursure','s','mln',3,'db3'); %denoise the signal
[z zloc]=max(sden); %identify the maximum value
zhalf=find(sden>=z*.9,1,'first');
z=trigt(zhalf)*1e9;   %report correction for time 0 value in nanoseconds
toff=47.9-z;    %time correction offset
time=((time.*10^9)+toff)';
%  clear file file2 zloc raw2 sden sden2 trigs z hm trigt;


%filter frequencies;
   %[a1 a2]=dwt(amp,'haar');
   af=fft(amp);
   af(1:3)=0;af(length(af)-3:length(af))=0;
   amp2=real(ifft(af)); clear af;
   
%Interpolate data with spline for better peak finding
   t2=(linspace(time(1),time(length(time)),4.*length(amp2)))';
   amp2=spline(time,amp,t2);
   time=t2;
   xd = wden(amp2,'heursure','s','mln',7,'dmey');
   amp2=xd;
%Find peaks

[high,low]=peakdet_v4(amp2,s,time);
%Determine displacement and velocity of peaks;
vh=0.775./(diff(high(:,1)));
vl=0.775./(diff(low(:,1)));
th=high(:,1);th(length(th))=[];
tl=low(:,1);tl(length(tl))=[];
alex=sortrows([th,vh;tl,vl]);

out_t=alex(:,1);
out_speed=alex(:,2);
out_time=time;
out_intensity=amp2;
end

