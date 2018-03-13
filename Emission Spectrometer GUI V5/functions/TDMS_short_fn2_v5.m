function [binned_time binned_spec_radiance  ] = TDMS_short_fn_v5( file,original_path,new_path,delay,unit_conversion,emission_filter_flag,calibration_filter_flag,starting_decade,ending_decade,binning_resolution,supressed_channels )
%TDMS_short_fn Shortens large TDMS files to small files and saves them into
%a new folder for accessing later

%Created shortned file name
name_short=[file(1:end-5),'_short.txt'];
full_file_name_short=fullfile(new_path,name_short);

raw_file=convertTDMS(0,fullfile(original_path,file));

A=extractfield(raw_file.Data.MeasuredData,'Data');

B=reshape(A,length(A)/33,33);

%Isolate time and flip voltage
time=B(:,5);B(:,5)=[];
PMT_out=B*-1;

%Baseline correction from the data in negative time
    for w=1:32
        baseline=mean(PMT_out([1:length(time)*.04],w));
        PMT_baseline_corrected(:,w)=PMT_out(:,w)-baseline;
    end

%     %Apply filters
% switch emiss_filter_flag
%     case 0
%         emiss_filter=0;
%     case 1
%         emiss_filter=0;
%     case 2
%         emiss_filter=0;
%     case 3
%         emiss_filter=0;
%     case 4
%         emiss_filter=0;
% end
% 
% switch cal_filter_flag
%     case 0
%         cal_filter=0;
%     case 1
%         cal_filter=0;
%     case 2
%         cal_filter=0;
%     case 3
%         cal_filter=0;
%     case 4
%         cal_filter=0;
% end

filter=ones(1,32);
    
%Correct for the 1 data point offset of card 1
tempcell=mat2cell(PMT_baseline_corrected,length(time),[4 28]);
tempcell{1,1}=[zeros(1,4);tempcell{1,1}];
tempcell{1,1}(end,:)=[];
PMT_time_corrected=[tempcell{1,1} tempcell{1,2}];

%Apply delays
time=time-delay;

%Remove negative data
vector=PMT_time_corrected(time>0,:);
time_vector=time(time>0);

%Set binning resolution and which decads to bin
res=binning_resolution;
str_dec=starting_decade;
end_dec=ending_decade;

%Begin data binning at 10 ns
binned_spec_radiance=vector(find(time_vector>0,1,'first'):find(time_vector<1*10^str_dec,1,'last'),:);
binned_time=time_vector(find(time_vector>0,1,'first'):find(time_vector<1*10^str_dec,1,'last'));

%Generate time range for each iteration
w=0;
for dec=str_dec:end_dec
timerange(1+w*res:(w+1)*res)=logspace(dec,dec+1,res);
w=w+1;
end

%Kill time range beyond the limits of our data
timerange=timerange(timerange<max(time_vector));


%Remove duplicates from the logtime generation
switch abs(str_dec)-abs(end_dec);
    case 6
timerange(res)=[];timerange(2*res-1)=[];timerange(3*res-2)=[];timerange(4*res-3)=[];timerange(5*res-4)=[];timerange(6*res-5)=[];
    case 5
timerange(res)=[];timerange(2*res-1)=[];timerange(3*res-2)=[];timerange(4*res-3)=[];timerange(5*res-4)=[];
    case 4
timerange(res)=[];timerange(2*res-1)=[];timerange(3*res-2)=[];timerange(4*res-3)=[];
    case 3
timerange(res)=[];timerange(2*res-1)=[];timerange(3*res-2)=[];
    case 2
timerange(res)=[];timerange(2*res-1)=[];
end

    %Generate vector of the indexes for binning ranges
    for logidx=1:length(timerange);
        index_vector(logidx)=find(time_vector>=timerange(logidx),1,'first');
    end

    %Generate temporary vector and build on output vectors
    for q=1:length(index_vector)-1;
        tempvector=vector(index_vector(q):index_vector(q+1),:);
        binned_spec_radiance(end+1,:)=mean(tempvector);
        temptime=time_vector(index_vector(q):index_vector(q+1));
        binned_time(end+1)=mean(temptime);
    end
    
    %Radiance calibration, filter, and voltage to amps conversions
    for q=1:size(binned_spec_radiance,1)
        binned_spec_radiance(q,:)=binned_spec_radiance(q,:).*unit_conversion./filter./50;
    end
    
    if supressed_channels ~=0;
%Supress relevant channels
binned_spec_radiance(:,supressed_channels)=1e-9;
    end
    
        %Remove very small data
    if binned_time(1,1)<1e-9
        binned_time(1)=[];
        binned_spec_radiance(1,:)=[];
    end
    
        %Condense output vectors
    output(:,1)=binned_time;
    output(:,[2:33])=binned_spec_radiance;
        

    
    %write text file with output vectors
    dlmwrite(full_file_name_short,output,' ');
    
end

