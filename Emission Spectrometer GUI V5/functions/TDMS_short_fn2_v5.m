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
switch emission_filter_flag
    case 0
        emiss_filter=ones(1,32);
    case 1
        emiss_filter=[2.22021012294715,2.18840061487304,2.28050090050574,2.22395699950986,2.15538712704376,2.17021904844884,2.15019386408174,2.18911043393480,2.14042885319131,2.12092026247839,2.10918767607242,2.05730096924894,2.06945514332392,2.02640359732492,2.01580500140237,2.03075017366455,2.00563520311422,2.00456730605755,1.99862587459952,1.97628167445285,1.97930107327330,1.95634363504716,1.96057644821198,1.96077846578524,1.91413206203148,1.89262982915889,0.867571561655801,1.86899219775357,1.85527558701627,1.84613287952972,1.84377084133934,2.20932386560421];
    case 2
        emiss_filter=ones(1,32);
    case 3
        emiss_filter=[19.3907025291586,19.5746394690905,19.4322103898503,19.0858027665230,18.0816751115418,18.1188506810188,17.7197928334547,18.4072453926705,17.6654800534518,17.2934718746698,17.5150530434812,16.7878855840659,16.6566057762737,16.4217983316711,15.9730755791696,16.0361562144140,15.6587658807397,15.3288760608162,15.3953125561919,15.1716634246300,14.8707637078099,15.1013134374759,14.6215955775044,14.6326512894228,14.2579711975102,14.1197568758042,-0.827617451532407,13.7316722747847,13.5147927580423,13.4003170609867,13.6069688191752,9.65288587328391];
    case 4
        emiss_filter=[2.52047401485964,2.43936120427156,2.38468834932386,2.32716071355890,2.18358010305433,2.15504169956444,2.08639427682733,2.11326108663841,2.05072896719365,1.99050743127860,1.97076455691803,1.92207078265933,1.91892157298665,1.85285211286551,1.88642134291774,1.92584694686741,1.90818714631405,1.91605234898119,1.92590472003531,1.91051882588537,1.92249418441883,1.89089848476122,1.89112787380383,1.89593493667419,1.84644973155697,1.81065950682353,0.954960219331430,1.74264512589156,1.69486455011534,1.65871012678058,1.57036544830753,1.67376168021530];
end

switch calibration_filter_flag
    case 0
        cal_filter=ones(1,32);
    case 1
        cal_filter=[2.22021012294715,2.18840061487304,2.28050090050574,2.22395699950986,2.15538712704376,2.17021904844884,2.15019386408174,2.18911043393480,2.14042885319131,2.12092026247839,2.10918767607242,2.05730096924894,2.06945514332392,2.02640359732492,2.01580500140237,2.03075017366455,2.00563520311422,2.00456730605755,1.99862587459952,1.97628167445285,1.97930107327330,1.95634363504716,1.96057644821198,1.96077846578524,1.91413206203148,1.89262982915889,0.867571561655801,1.86899219775357,1.85527558701627,1.84613287952972,1.84377084133934,2.20932386560421];
    case 2
        cal_filter=ones(1,32);
    case 3
        cal_filter=[19.3907025291586,19.5746394690905,19.4322103898503,19.0858027665230,18.0816751115418,18.1188506810188,17.7197928334547,18.4072453926705,17.6654800534518,17.2934718746698,17.5150530434812,16.7878855840659,16.6566057762737,16.4217983316711,15.9730755791696,16.0361562144140,15.6587658807397,15.3288760608162,15.3953125561919,15.1716634246300,14.8707637078099,15.1013134374759,14.6215955775044,14.6326512894228,14.2579711975102,14.1197568758042,-0.827617451532407,13.7316722747847,13.5147927580423,13.4003170609867,13.6069688191752,9.65288587328391];
    case 4
        cal_filter=[2.52047401485964,2.43936120427156,2.38468834932386,2.32716071355890,2.18358010305433,2.15504169956444,2.08639427682733,2.11326108663841,2.05072896719365,1.99050743127860,1.97076455691803,1.92207078265933,1.91892157298665,1.85285211286551,1.88642134291774,1.92584694686741,1.90818714631405,1.91605234898119,1.92590472003531,1.91051882588537,1.92249418441883,1.89089848476122,1.89112787380383,1.89593493667419,1.84644973155697,1.81065950682353,0.954960219331430,1.74264512589156,1.69486455011534,1.65871012678058,1.57036544830753,1.67376168021530];
end

filter=cal_filter/emiss_filter;
    
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

