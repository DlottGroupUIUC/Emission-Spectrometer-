function [ delay ] = peakdelay( file,original_path,target,threshold)
%PeakDelay10ns function is designed to extract the time point where the
%nanosecond peak in standard molecular explosive emission is at 10 ns.
%  
%import data
raw_file=convertTDMS(0,fullfile(original_path,file));

time=raw_file.Data.MeasuredData(8).Data;
line_out=-raw_file.Data.MeasuredData(25).Data;

%Clip to two microseconds
line_out(time<-2e-7)=[];
time(time<-2e-7)=[];
line_out(time>2e-6)=[];
time(time>2e-6)=[];

rms=sqrt(mean(line_out(1:250).^2));

[~,locs] = findpeaks(line_out,'MinPeakHeight',threshold*rms,...
'MinPeakDistance',50e-9);

[~, loc_max]=max(line_out(locs(1:5)));

delay=time(locs(loc_max))-target;
end

