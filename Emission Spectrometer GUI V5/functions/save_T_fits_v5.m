function [] = save_with_flags_v5( path,data_block,saved_file_name)
% Temperature_vector_save Saves the output temperature vector to an ascii
%file with proper headers such that it can be imported into origin simply

%Generate file name
file_name=fullfile(path,saved_file_name);

%Replace empty cells with NaN
empty_cells=cellfun('isempty',data_block);
data_block(empty_cells)={NaN};

%Headers
for idx=0:size(data_block,1)-2
    
    %set header for wavelength once
    if idx==0
        
        %Title
        hdr{1}='wavelength';
        hdr{2}='spectral radiance';
        hdr{3}='model';
        hdr{4}='fit parameters';
        
        %Units
        hdr2{1}='nm';
        hdr2{2}='W/sr m2 nm';
        hdr2{3}='W/sr m2 nm';
        hdr2{4}='T(K) \phi(unitless)';
        
        %Comments
        hdr3{1}='';
        hdr3{2}=data_block{2,4};
        hdr3{3}='';
        hdr3{4}='';
        
        %wavelength data
        data(:,4*idx+1)=data_block{1,1};
                
        %spectral radiance data
        data(:,4*idx+2)=data_block{2,1};
        
        %model data
        data(:,4*idx+3)=data_block{2,2};
        
        %model parameters data
        data=padconcatenation_v4(data,data_block{idx+2,3}',2);
         
    else 
       
        %Title
        hdr{3*idx+2}='spectral radiance';
        hdr{3*idx+3}='model';
        hdr{3*idx+4}='fit parameters';
        
        %Units
        hdr2{3*idx+2}='W/sr m2 nm';
        hdr2{3*idx+3}='W/sr m2 nm';
        hdr2{3*idx+4}='T(K) \phi(unitless)';
        
        %Comments
        hdr3{3*idx+2}=data_block{idx+2,4};
        hdr3{3*idx+3}='';
        hdr3{3*idx+4}='';
                
        %spectral radiance data
        data=padconcatenation_v4(data,data_block{idx+2,1}',2);
        
        %model data
        data=padconcatenation_v4(data,data_block{idx+2,2}',2);
        
        %model parameters data
        data=padconcatenation_v4(data,data_block{idx+2,3}',2);
        
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




