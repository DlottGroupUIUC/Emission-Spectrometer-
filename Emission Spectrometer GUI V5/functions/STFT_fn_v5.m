function [out_t,out_speed,out_time,out_intensity,out_SG_v,out_SG_t,out_SG_spectrogram] = STFT_fn_v5(file_name,path,channel,PDV_y_cutoff)
%pdverror_fn A functional form of the PDVerror script for use with
%GUIs

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

%Load Trigger Data on Channel 3
file2=importdata(fullfile(path,[file_name,'_ch3.txt']));
trigt=time;  %separate data file into time
trigs=file2.data(:,2);  %and signal variables
[z zloc]=max(trigs); %identify the maximum value
zhalf=find(trigs>=z*.9,1,'first');
z=trigt(zhalf)*1e9;   %report correction for time 0 value in nanoseconds
toff=47.9-z;    %time correction offset

vcutoff=PDV_y_cutoff; %Lower velocity cutoff

%Determine frequency of spectrogram
fs=1./(time(2)); %sampling frequency
g=160;
r=floor(fs./(g.*10^6)); %sets r as the number of points needed for given resolution
clear g;
   af=fft(amp);          %can replace with using a1 for amplification of signal and faster processing but lower velocity range
   amp2=real(ifft(af));clear af;           
%Set amplitude graphing bounds
   ampmax=round((max(abs(min(amp2)),max(amp2))+0.05)*10)/10;ampmin=-ampmax; %makes images cleaner when exporting
%Calculate Spectrogram using Hanning Window
[R,f,t]=spectrogram(amp2,hann(floor(r)),floor(r)-1,5.*floor(r),fs);
v=f.*0.775./1e9;
clear f;
t=((t.*10^9)+toff)';
R=abs(R);
%Convert time from seconds to nanoseconds
time=time.*10^9+toff;
%Apply low frequency filter
filter=length(v(v<vcutoff));
R(1:filter,:)=0;
clear filter;
clear vcutoff;

% %Plot spectrogram in secondary window
% figure(2);imagesc(time,v,R);set(gca,'ydir','normal'); ylabel('velocity (km/s)'); ylim([0 2]);
% xlim([0, time(length(time))]); xlabel('time (ns)');

%Create single line spectrogram
[mx locs]=max(R,[],1);
velocities=v(locs);

t2=(linspace(time(1),time(length(time)),length(amp2)))';

out_t=t;
out_speed=velocities;
out_time=t2;
out_intensity=amp2;

out_SG_v=v;
out_SG_t=t;
out_SG_spectrogram=R;


end

