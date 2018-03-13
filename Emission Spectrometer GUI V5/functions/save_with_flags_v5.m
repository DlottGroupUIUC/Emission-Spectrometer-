function [] = save_with_flags_v5( path,file_list,data_block,saved_file_name,save_type_flag )
% Temperature_vector_save Saves the output temperature vector to an ascii
%file with proper headers such that it can be imported into origin simply

%Generate file name
file_name=fullfile(path,saved_file_name);

%Replace empty cells with NaN
empty_cells=cellfun('isempty',data_block);
data_block(empty_cells)={NaN};

%Headers
for idx=0:length(file_list)-1
    
    switch save_type_flag
        
        % Case 1 refers to saving impact data and has a different format
        case 1
            
            if idx==0
                %set number of items to save
                num=4;
                
                %Long name
                hdr{1}='shot number';
                hdr{2}='speed';
                hdr{3}='speed error';
                hdr{4}='impact timing';
                
                %Units
                hdr2{1}='';
                hdr2{2}='km/s';
                hdr2{3}='km/s';
                hdr2{4}='s';
                
                %Comments
                hdr3{1}='';
                hdr3{2}='';
                hdr3{3}='';
                hdr3{4}='';
                
                %Unwrap data cell array
                %grab shot number
                data(:,1)=cellfun(@str2num,data_block(:,1));
                %speed
                data(:,2)=data_block{:,2};
                %speed error
                data(:,3)=data_block{:,3};
                %impact time
                data(:,4)=data_block{:,4};
                
            end
            
        % case 2 refers to saving radiance. Its format is consistent with
        %the next several cases
        case 2
            
            %set number of items to save
            num=2;
            
            %Long name
            hdr{num*idx+1}='time';
            hdr{num*idx+2}='radiance';
            
            %Units
            hdr2{num*idx+1}='s';
            hdr2{num*idx+2}='W/sr m2';
            
            %Comments
            hdr3{num*idx+1}=file_list{idx+1};
            hdr3{num*idx+2}=file_list{idx+1};
            
            %Unwrap data cell array
            if idx==0
                
                %time data
                data(:,4*idx+1)=data_block{1,2};
                
                %radiance data
                data(:,4*idx+2)=data_block{1,3};
                
                
            else
                
                %time data
                data=padconcatenation_v4(data,data_block{idx+1,2},2);
                
                %radiance data
                data=padconcatenation_v4(data,data_block{idx+1,3},2);
                
                
            end
            
        % case 3 deals with saving Z fits.
        case 3
            
            %set number of items to save
            num=3;
            
            %Long name
            hdr{num*idx+1}='time';
            hdr{num*idx+2}='temperature';
            hdr{num*idx+3}='error';
            
            %Units
            hdr2{num*idx+1}='s';
            hdr2{num*idx+2}='K';
            hdr2{num*idx+3}='K';
            
            
            %Comments
            hdr3{num*idx+1}=file_list{idx+1};
            hdr3{num*idx+2}=file_list{idx+1};
            hdr3{num*idx+3}=file_list{idx+1};
            
            %Unwrap data cell array
            if idx==0
                
                %time data
                data(:,4*idx+1)=data_block{1,2};
                
                %temperature data
                data(:,4*idx+2)=data_block{1,5};
                
                %temperature error data
                data(:,4*idx+3)=data_block{1,6};
                
                
            else
                
                %time data
                data=padconcatenation_v4(data,data_block{idx+1,2},2);
                
                %temperature data
                data=padconcatenation_v4(data,data_block{idx+1,5},2);
                
                %temperature error data
                data=padconcatenation_v4(data,data_block{idx+1,6},2);
                
                
            end
            
        % case 4 deals with saving gray temperature fits.
        case 4
            
            %set number of items to save
            num=3;
            
            %Long name
            hdr{num*idx+1}='time';
            hdr{num*idx+2}='temperature';
            hdr{num*idx+3}='error';
            
            %Units
            hdr2{num*idx+1}='s';
            hdr2{num*idx+2}='K';
            hdr2{num*idx+3}='K';
            
            
            %Comments
            hdr3{num*idx+1}=file_list{idx+1};
            hdr3{num*idx+2}=file_list{idx+1};
            hdr3{num*idx+3}=file_list{idx+1};
            
            %Unwrap data cell array
            if idx==0
                
                %time data
                data(:,4*idx+1)=data_block{1,2};
                
                %temperature data
                data(:,4*idx+2)=data_block{1,7};
                
                %temperature error data
                data(:,4*idx+3)=data_block{1,8};
                
                
            else
                
                %time data
                data=padconcatenation_v4(data,data_block{idx+1,2},2);
                
                %temperature data
                data=padconcatenation_v4(data,data_block{idx+1,7},2);
                
                %temperature error data
                data=padconcatenation_v4(data,data_block{idx+1,8},2);
                
                
            end
            
        % case 5 deals with saving gray phi fits.
        case 5
            
            %set number of items to save
            num=3;
            
            %Long name
            hdr{num*idx+1}='time';
            hdr{num*idx+2}='phi';
            hdr{num*idx+3}='error';
            
            %Units
            hdr2{num*idx+1}='s';
            hdr2{num*idx+2}=' ';
            hdr2{num*idx+3}=' ';
            
            
            %Comments
            hdr3{num*idx+1}=file_list{idx+1};
            hdr3{num*idx+2}=file_list{idx+1};
            hdr3{num*idx+3}=file_list{idx+1};
            
            %Unwrap data cell array
            if idx==0
                
                %time data
                data(:,4*idx+1)=data_block{1,2};
                
                %temperature data
                data(:,4*idx+2)=data_block{1,9};
                
                %temperature error data
                data(:,4*idx+3)=data_block{1,10};
                
                
            else
                
                %time data
                data=padconcatenation_v4(data,data_block{idx+1,2},2);
                
                %temperature data
                data=padconcatenation_v4(data,data_block{idx+1,9},2);
                
                %temperature error data
                data=padconcatenation_v4(data,data_block{idx+1,10},2);
                
                
            end
            
            
    end
    
end

% the engine
txt=sprintf('%s\t',hdr{:});
txt(end)='';
txt2=sprintf('%s\t',hdr2{:});
txt2(end)='';
txt3=sprintf('%s\t',hdr3{:});
txt3(end)='';
dlmwrite(file_name,txt,'');
dlmwrite(file_name,txt2,'-append','delimiter','');
dlmwrite(file_name,txt3,'-append','delimiter','');
dlmwrite(file_name,data,'-append','delimiter','\t');




