function meta_data = snn_dispatch_meta_data( file_name, section_name, num_fields )
% snn_dispatch_meta_data: reads meta data from m-file
%
% meta_data = snn_dispatch_meta_data( file_name, section_name, num_fields )
%
% Reads meta data from m-file and convers each numerical
% field to numerical values.
%
% see also
%   <a href="matlab:help snn_parse_meta_data">snn_parse_meta_data</a>
%   

    meta_data = snn_parse_meta_data( file_name, section_name, num_fields );

    for p = 1:(num_fields+1):length(meta_data)        
        for i=1:num_fields
            
            d_value = meta_data{p+i};
            
            if ( d_value(1) == '''' ) && ( d_value(end) == '''' )
                meta_data{p+i} = d_value(2:end-1);
            else
                num_value = str2num( d_value );

                if ~isempty( num_value )
                    meta_data{p+i} = num_value;
                end
            end
        end
    end
end
