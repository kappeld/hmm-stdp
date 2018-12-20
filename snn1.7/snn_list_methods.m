function snn_list_methods( type )
% snn_list_methods: lists all SNN methods.
%
% snn_list_methods( type )
% snn_list_methods
%
% Displays a list of all available SNN methods of specified
% type and a brief text descibing them. If no type is given
% all three types are displayed. Directories that contain
% SNN methods must be added to the SNN search directories
% using the <a href="matlab:help snn_include">snn_include</a> function.
%
% input
%   type: Can be 'train', 'sample' or 'performance'.
%         If nothing is passed all types are displayed.
%
% see also
%   <a href="matlab:help snn_include">snn_include</a>
%
% David Kappel 17.11.2010
%

    if (nargin==0)        
        snn_list_methods( 'train' );
        snn_list_methods( 'sample' );
        snn_list_methods( 'performance' );
    else
        if ~ischar( type )
            error( 'Unexpected input type!' );
        end
        
        if ~strcmp( type, 'train' ) && ...
           ~strcmp( type, 'sample' ) && ...
           ~strcmp( type, 'performance' )
            error( 'Unknown method type ''%s''!', type );
        end

        search_dirs = snn_options( 'search_dirs' );

        fprintf( '%s methods:\n', type );
        
        for p = 1:length(search_dirs)
            
            path_info = what( search_dirs{p} );
        
            if isempty( path_info ) || isempty( path_info.path )
                continue;
            end
        
            method_list = dir( [ path_info.path, '/snn_', type, '_*.m' ] );
        
            for i=1:length( method_list )

                def_str = method_list(i).name( (6+length(type)):end-2 );

                m_lines = textread( method_list(i).name, '%s', 'delimiter','\n' );

                description = '';

                for line_no=1:length(m_lines)

                    if ( m_lines{line_no}(1) == '%' )
                        description = m_lines{line_no}( 2:end );
                        break;
                    end
                end

                fprintf( '  ''<a href = "matlab:help %s">%s</a>''%s%s\n', ...
                         method_list(i).name, def_str, repmat(' ', 20-length(def_str), 1 ), ...
                         description );
            end
        end
        
        fprintf( '\n' );
    end
end
