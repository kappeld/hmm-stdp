function [tokens,start_pos,end_pos] = snn_parse_args( in_str )
% Dispatch an argument string
%
% tokens = snn_parse_args( in_str )
%
% Parse an rgument of the form [<arg1>,<arg2>,...,<argN>] into
% a cell array of string tokens <arg1>,...,<argN>.
%
% 17.11.2010
%

    in_str = in_str( in_str ~= ' ' );
    start_pos = strfind( in_str, '[' );
    end_pos = strfind( in_str, ']' );
    
    if ( length( start_pos ) ~= length( end_pos ) )
        error( 'Unbalanced braces pair in allocator comamnd ''%s''!', in_str );
    end
    
    if isempty( start_pos ) || isempty( end_pos )
        
        start_pos = length(in_str)+1;
        end_pos = length(in_str)+1;
        tokens = {};
        return;
    end

    start_pos = start_pos(1);
    end_pos = end_pos(end);
    
    tokens = strread( in_str( start_pos+1:end_pos-1 ), '%s', 'delimiter', ',' );
end
