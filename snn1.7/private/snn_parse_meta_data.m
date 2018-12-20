function params = snn_parse_meta_data( file_name, section_name, num_fields )
% SNN_PARSE_META_DATA reads the meta data from an m-file
%
% params = snn_parse_meta_data( file_name, section_name, num_fields )
%
% Reads meta data from m-file. The meta data within the m-files
% documentation text is used to adjust the the network structure. Meta data
% is organised in sections. A meta data section starts with a line:
% @<section name>: and ends with an empty comment line (a line that only
% contains a '%'). The meta data sections may contrain any number of
% parameter value pairs separated with white spaces and may be followed by
% some describing text. Eeach line contrains a single parameter followed by
% num_fields values. Meta data is used to adjust the behaviour of the
% network. Each SNN method should define its default parameters in the
% 'default parameters' meta data section.
%
% example:
%
% @default parameters:
%   eta           0.5              Here I define the default learn rate
%   train_method  st               Here I define the default sample method
%   some_matrix   [ 0, 0, 10, 10 ] Here I define some matrix
%   some_string   'some text'      Here I define some string
%

    meta_data_lines = textread( file_name, '%s', 'delimiter','\n' );
        
    sections = strmatch( [ '% @', section_name, ':' ], meta_data_lines );
    
    params = cell( 0, num_fields+1 );
    
    for sec = 1:length(sections)
        
        cur_line = sections(sec);
        
        while true
           
            cur_line = cur_line+1;
            
            toks = strread( meta_data_lines{cur_line}(2:end), '%s' );
            
            if ( isempty(toks) )
                break;
            end
            
            if ( toks{1}(1) == '@')
                break;
            end
            
            pars = cell( 1, num_fields+1 );
            par_depth = 0;
            str_tok = false;
            idx = 1;
            
            for i=1:length(toks)
                
                if ( toks{i}(1) == '''' )
                    str_tok = true;
                end

                pars{idx} = [ pars{idx}, toks{i} ];
                
                if str_tok
                    if ( toks{i}(end) == '''' )
                        str_tok = false;
                    else
                        pars{idx} = [ pars{idx}, ' ' ];                        
                    end
                else                
                    par_depth = par_depth + length( strfind( toks{i}, '[' ) ) - ...
                                            length( strfind( toks{i}, ']' ) );
                end

                if ( par_depth == 0 ) && ~str_tok
                    
                    idx = idx+1;
                    
                    if ( idx > num_fields+1 )
                        break;
                    end
                end
            end
            
            if ( par_depth ~= 0 ) || str_tok
                error( [ 'Error while parsing meta data in file ''%s'' at line %u!\n', ...
                         'Unbalanced braces or quotes!' ], file_name, cur_line );
            end

            params = { params{:}, pars{:} };
            
        end
    end
end
