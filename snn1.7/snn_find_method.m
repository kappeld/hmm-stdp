function method_name = snn_find_method( type, pattern )
% snn_find_method: returns the full name of spcified SNN method
%
% method_name = snn_find_method( type, pattern )
% method_name = snn_find_method( type )
%
% Retruns the full name of the specified SNN method.
% Searches all snn search directories and returns the first
% method name that matches the pattern. If no method
% was found an error is thrown.
%
% input
%   type:        The snn method type. Can be 'train',
%                'sample' or 'performance'.
%   pattern:     A name pattern string that may contain
%                place holders like '*'. Default is '*'.
%
% output
%   method_name: The full method name.
%
% see also
%   <a href="matlab:help snn_include">snn_include</a>
%
% Daid Kappel 17.11.2010
%

    if (nargin<1)
        errro( 'Not enough input arguments!' );
    end
    
    if (nargin<2)
        method_name = '*';
    end
    
    if ~ischar( type )
        error( 'Unexpected input type for argument ''type''!' );
    end

    if ~ischar( pattern )
        error( 'Unexpected input type for argument ''method_name''!' );
    end

    
    if ~strcmp( type, 'train' ) && ...
       ~strcmp( type, 'sample' ) && ...
       ~strcmp( type, 'performance' )
        error( 'Unexpected input ''%s''.', type );
    end

    search_dirs = snn_options( 'search_dirs' );

    method = [ 'snn_', type, '_', pattern ];

    for p = 1:length(search_dirs)

        path_info = what( search_dirs{p} );

        if isempty( path_info ) || isempty( path_info.path )
            continue;
        end

        method_list = dir( [ path_info.path, '/', method, '.m' ] );
        
        if ~isempty( method_list )
            method_name = method_list(1).name;
            method_name = method_name(1:end-2);
            return;
        end
    end

    error( 'No %s method found matching pattern ''%s''!', type, pattern );
end
