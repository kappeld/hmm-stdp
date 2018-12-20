function snn_print( str, varargin )
%  snn_print( str, varargin )
%
% David Kappel
% 01.12.2011
%

    global snn_options__;
    
    if isempty( snn_options__ )        
        snn_options__ = struct;        
    end

    if ( nargin > 0 )
        
        if ~ischar( op_name )
            error( 'Unexpected argument type for ''op_name''!' );
        end
        
        if isfield( snn_options__, op_name )
            value = snn_options__.(op_name);
        else
            value = [];
        end
        
        if ( nargin > 1 )            
            if isempty( new_value )
                if isfield( snn_options__, op_name )
                    snn_options__ = rmfield( snn_options__, op_name );
                end
            else
                snn_options__.(op_name) = new_value;
            end
        end
    else        
        value = snn_options__;
    end
end
